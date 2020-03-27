import numpy as np
import pandas as pd
from ibl_pipeline.utils import psychofit
import query_scan_mice as query
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil
from scipy import stats
import matplotlib as mpl
import functools
import numpy as np
import scipy.optimize
from scipy.special import erf, erfc

def mle_fit_psycho(data, P_model='weibull', parstart=None, parmin=None, parmax=None, nfits=50):
    """
    Maximumum likelihood fit of psychometric function.

    Args:
        data: 3 x n matrix where first row corrsponds to stim levels (%), 
            the second to number of trials for each stim level (int),
            the third to proportion correct (float between 0 and 1)
        P_model: The psychometric function. Possibilities include 'weibull'
            (DEFAULT), 'weibull50', 'erf_psycho' and 'erf_psycho_2gammas'
        parstart: Non-zero starting parameters, used to try to avoid local
            minima.  The parameters are [threshold, slope, gamma], or if
            using the 'erf_psycho_2gammas' model append a second gamma value.
            Recommended to use a value > 1.
            If None, some reasonable defaults are used.
        parmin: Minimum parameter values.  If None, some reasonable defaults 
            are used
        parmax: Maximum parameter values.  If None, some reasonable defaults 
            are used
        nfits: the number of fits

    Returns:
        pars: The parameters from the best of the fits
        L: The likliehood of the best fit
        
    Raises:
        TypeError: data must be a list or numpy array
        ValueError: data must be m by 3 matrix

    Examples:
        Below we fit a Weibull function to some data:
        
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> cc = np.array([-8., -6., -4., -2.,  0.,  2.,  4.,  6.,  8.]) # contrasts
        >>> nn = np.full((9,),10) # number of trials at each contrast
        >>> pp = np.array([5., 8., 20., 41., 54., 59., 79., 92., 96])/100 # proportion "rightward"
        >>> pars, L = mle_fit_psycho(np.vstack((cc,nn,pp)), 'erf_psycho')
        >>> plt.plot(cc, pp, 'bo', mfc='b')
        >>> plt.plot(np.arange(-8,8,0.1), erf_psycho(pars,np.arange(-8,8,0.1)), '-b')
        
    Information:
        1999-11 FH wrote it
        2000-01 MC cleaned it up
        2000-04 MC took care of the 50% case
        2009-12 MC replaced fmins with fminsearch
        2010-02 MC, AZ added nfits
        2013-02 MC+MD fixed bug with dealing with NaNs
        2018-08 MW ported to Python
    """
    # Input validation
    if isinstance(data, (list, tuple)):
        data = np.array(data)
    elif not isinstance(data, np.ndarray):
        raise TypeError('data must be a list or numpy array')

    if data.shape[0] != 3:
        raise ValueError('data must be m by 3 matrix')

    if parstart is None:
        parstart = np.array([np.mean(data[0,:]), 3., .05])
    if parmin is None:
        parmin = np.array([np.min(data[0,:]), 0., 0.])
    if parmax is None:
        parmax = np.array([np.max(data[0,:]), 10., .4])

    # find the good values in pp (conditions that were effectively run)
    ii = np.isfinite(data[2,:]);

    likelihoods = np.zeros(nfits,);
    pars = np.empty((nfits,parstart.size))
    
    f = functools.partial(neg_likelihood, data=data[:,ii], 
                          P_model=P_model, parmin=parmin, parmax=parmax)
    for ifit in range(nfits):
        pars[ifit,:] = scipy.optimize.fmin(f, parstart, disp=False)
        parstart = parmin + np.random.rand(parmin.size) * (parmax-parmin)
        likelihoods[ifit] = - neg_likelihood(pars[ifit,:], data[:,ii], P_model, parmin, parmax)

    # the values to be output
    L = likelihoods.max()
    iBestFit = likelihoods.argmax()
    return pars[iBestFit,:], L

def neg_likelihood(pars, data, P_model='psych_cumNorm', parmin=None, parmax=None):
    """
    Negative likelihood of a psychometric function.

    Args:
        pars: Model parameters [threshold, slope, gamma], or if
            using the 'erf_psycho_2gammas' model append a second gamma value.
        data: 3 x n matrix where first row corrsponds to stim levels (%), 
            the second to number of trials for each stim level (int),
            the third to proportion correct (float between 0 and 1)
        P_model: The psychometric function. Possibilities include 'weibull'
            (DEFAULT), 'weibull50', 'erf_psycho' and 'erf_psycho_2gammas'
        parmin: Minimum bound for parameters.  If None, some reasonable defaults 
            are used
        parmax: Maximum bound for parameters.  If None, some reasonable defaults 
            are used

    Returns:
        l: The likliehood of the parameters.  The equation is:
            - sum(nn.*(pp.*log10(P_model)+(1-pp).*log10(1-P_model)))
            See the the appendix of Watson, A.B. (1979). Probability
            summation over time. Vision Res 19, 515-522.
        
    Raises: 
        ValueError: invalid model, options are "weibull", 
                    "weibull50", "erf_psycho" and "erf_psycho_2gammas"
        TypeError: data must be a list or numpy array
        ValueError data must be m by 3 matrix
        
    Information:
        1999-11 FH wrote it
        2000-01 MC cleaned it up
        2000-07 MC made it indep of Weibull and added parmin and parmax
        2018-08 MW ported to Python
    """
    # Validate input
    if isinstance(data, (list, tuple)):
        data = np.array(data)
    elif not isinstance(data, np.ndarray):
        raise TypeError('data must be a list or numpy array')
        
    if parmin is None:
        parmin = np.array([.005, 0., 0.])
    if parmax is None:
        parmax = np.array([.5, 10., .25])
        
    if data.shape[0] == 3:
        xx = data[0,:]
        nn = data[1,:]
        pp = data[2,:]
    else:
        raise ValueError('data must be m by 3 matrix')

    # here is where you effectively put the constraints.
    if (any(pars < parmin)) or (any( pars > parmax)):
        l = 10000000
        return l

    dispatcher={'psych_cumNorm': psych_cumNorm}
    try:
        probs = dispatcher[P_model](pars,xx)
    except KeyError:
        raise ValueError('invalid model, options are "weibull", '+
                         '"weibull50", "erf_psycho" and "erf_psycho_2gammas"')

    assert (max(probs)<=1) or (min(probs) >= 0),'At least one of the probabilities is not between 0 and 1'

    probs[probs==0]=np.finfo(float).eps
    probs[probs==1]=1-np.finfo(float).eps

    l = - sum(nn*(pp*np.log(probs)+(1-pp)*np.log(1-probs)))
    return l

def psych_cumNorm(params, x):

    alpha = params[0]
    beta = params[1]
    gamma = params[2]
    lamb = params[3]
    y = gamma + (1 - gamma - lamb) * .5 * erfc(-beta *(x-alpha) / np.sqrt(2))
    return y


bregma = [228.5, 190]
pixelsize = .025  # mm
allenOutline = np.load(r'F:\allen_dorsal_outline')
data = query.align_laser2behavior(['CSK-scan-014','CSK-scan-015','CSK-scan-016','CSK-scan-019'])

bigData = [pd.concat(data)]
data.append(bigData[0])
spots = data[0].groupby(['laserPosX', 'laserPosY']).size().reset_index().rename(columns={0: 'count'})
spotBias = [[[] for subject in range(len(data))] for spots in range(len(spots))]
spotSlope = [[[] for subject in range(len(data))] for spots in range(len(spots))]
spotLapseHigh = [[[] for subject in range(len(data))] for spots in range(len(spots))]
spotLapseLow = [[[] for subject in range(len(data))] for spots in range(len(spots))]
animalFits = []
subIdx = 0

plotContrasts = [1, .25,  .125, .0625, 0]  # 1, .25, .125, list of contrasts to use in spot plot

visLeftSpots = np.array([[-2.5,-1.5],[-3.5,-1.5],[-2.5,-2.5],[-3.5,-2.5]])
visRightSpots = np.array([[2.5,-1.5],[3.5,-1.5],[2.5,-2.5],[3.5,-2.5]])

moLeftSpots = np.array([[-1.5,.5],[-1.5,1.5],[-1.5,2.5],[-2.5,1.5]])
moRightSpots = np.array([[1.5,.5],[1.5,1.5],[1.5,2.5],[2.5,1.5]])

for subject in data:
    shuffledData = subject.copy()
    shuffledData['laserPosX'] = np.random.permutation(shuffledData.loc[:, 'laserPosX'])
    shuffledData['laserPosY'] = np.random.permutation(shuffledData.loc[:, 'laserPosY'])

    spots = subject.groupby(['laserPosX', 'laserPosY']).size().reset_index().rename(
        columns={0: 'count'})
    contrastSet = np.unique(subject['signed_contrast'])
    subject['CW'] = subject['trial_response_choice']=='CCW'
    psychoSpotData = [[contrastSet, [], []] for spot in range(len(spots))]
    spotData = [[[], [], [], [], []] for spot in range(len(spots))]
    pVals = [[[], [], [], [], [], []] for spot in range(len(spots))]
    controlData = [[], [], [], [], []]
    means = [[[], [], [], [], [], []] for spot in range(len(spots))]
    spotFits = []
    goLeft = 0
    goRight = 1
    noGo = 2
    RT = 3
    correct = 4

    # f, axes = plt.subplots(13, 10, figsize=(14, 14), sharex=True, sharey=True)
    # sns.despine(left=True, bottom=True)
    controlSpots = list(spots[spots['laserPosY'] <-6].index)
    for i in range(len(spots)):
        spot = spots.iloc[i, [0, 1]]
        tempX = subject[subject['laserPosX'] == spot[0]]
        tempSpot = tempX[tempX['laserPosY'] == spot[1]]

        for contrast in contrastSet:
            byContrast = tempSpot[tempSpot['signed_contrast'] == contrast]
            psychoSpotData[i][1].append(len(byContrast))
            # psychoSpotData[i][2].append(np.mean(byContrast['trial_feedback_type']))
            psychoSpotData[i][2].append(np.mean(byContrast['CW']))
            if  abs(contrast) in plotContrasts:
                spotData[i][goLeft].append(byContrast['trial_response_choice'] == 'CCW')  # go left
                spotData[i][goRight].append(byContrast['trial_response_choice'] == 'CW')  # goright
                spotData[i][noGo].append(byContrast['trial_response_choice'] == 'No Go')  # no go
                spotData[i][RT].append(byContrast['trial_response_time'] - byContrast[
                    'trial_go_cue_trigger_time'])
                byContrast[byContrast['trial_feedback_type'] == -1] = 0
                spotData[i][correct].append(byContrast['trial_feedback_type'])  # correct
        # Compiling lists
        spotData[i][goLeft] = [item for sublist in spotData[i][goLeft] for item in sublist]
        spotData[i][goRight] = [item for sublist in spotData[i][goRight] for item in sublist]
        spotData[i][noGo] = [item for sublist in spotData[i][noGo] for item in sublist]
        spotData[i][RT] = [item for sublist in spotData[i][RT] for item in sublist]
        spotData[i][correct] = [item for sublist in spotData[i][correct] for item in sublist]


    controlData[goLeft] = spotData[controlSpots[0]][goLeft] + spotData[controlSpots[1]][goLeft]
    controlData[goRight] = spotData[controlSpots[0]][goRight] + spotData[controlSpots[1]][goRight]
    controlData[noGo] = spotData[controlSpots[0]][noGo] + spotData[controlSpots[1]][noGo]
    controlData[RT] = spotData[controlSpots[0]][RT] + spotData[controlSpots[1]][RT]
    controlData[correct] = spotData[controlSpots[0]][correct] + spotData[controlSpots[1]][correct]

    ## Getting Pvals for each spot, and sorting out spots into visLeft and Right
    numTrials = 0
    for i in range(len(spots)):
        t, pVals[i][goLeft] = stats.ttest_ind(spotData[i][goLeft], controlData[goLeft])
        t, pVals[i][goRight] = stats.ttest_ind(spotData[i][goRight], controlData[goRight])
        t, pVals[i][noGo] = stats.ttest_ind(spotData[i][noGo], controlData[noGo])
        t, pVals[i][RT] = stats.ttest_ind(spotData[i][RT], controlData[RT])
        t, pVals[i][correct] = stats.ttest_ind(spotData[i][correct], controlData[correct])
        means[i][goLeft] = np.nanmean(spotData[i][goLeft])
        means[i][goRight] = np.nanmean(spotData[i][goRight])
        means[i][noGo] = np.nanmean(spotData[i][noGo])
        means[i][RT] = np.nanmedian(spotData[i][RT])
        means[i][correct] = np.nanmean(spotData[i][correct])
        numTrials+=len(spotData[i][goLeft])

    ## making psychometrics and other plots for vis inactivations vs control spots
    visLeftPsycho = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]
    visRightPsycho = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]

    moRightPsycho = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]
    moLeftPsycho = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]

    controlSpotPsycho = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]
    for i in range(len(spots)):
        spotX = spots.iloc[i, [0]][0]
        spotY = spots.iloc[i, [1]][0]
        for leftSpots in range(len(visLeftSpots)):
            if spotX == visLeftSpots[leftSpots][0] and spotY == visLeftSpots[leftSpots][1]:
                visLeftPsycho[0] = psychoSpotData[i][0]
                visLeftPsycho[1] = [i+j for i,j in zip(psychoSpotData[i][1],visLeftPsycho[1])]
                visLeftPsycho[2] = [np.nanmean([i,j]) for i,j in zip(psychoSpotData[i][2],visLeftPsycho[2])]
        for rightSpots in range(len(visRightSpots)):
            if spotX == visRightSpots[rightSpots][0] and spotY == visRightSpots[rightSpots][1]:
                visRightPsycho[0] = psychoSpotData[i][0]
                visRightPsycho[1] = [i+j for i,j in zip(psychoSpotData[i][1],visRightPsycho[1])]
                visRightPsycho[2] = [np.nanmean([i,j]) for i,j in zip(psychoSpotData[i][2],visRightPsycho[2])]

        for leftSpots in range(len(moLeftSpots)):
            if spotX == moLeftSpots[leftSpots][0] and spotY == moLeftSpots[leftSpots][1]:
                moLeftPsycho[0] = psychoSpotData[i][0]
                moLeftPsycho[1] = [i+j for i,j in zip(psychoSpotData[i][1],moLeftPsycho[1])]
                moLeftPsycho[2] = [np.nanmean([i,j]) for i,j in zip(psychoSpotData[i][2],moLeftPsycho[2])]
        for rightSpots in range(len(moRightSpots)):
            if spotX == moRightSpots[rightSpots][0] and spotY == moRightSpots[rightSpots][1]:
                moRightPsycho[0] = psychoSpotData[i][0]
                moRightPsycho[1] = [i+j for i,j in zip(psychoSpotData[i][1],moRightPsycho[1])]
                moRightPsycho[2] = [np.nanmean([i,j]) for i,j in zip(psychoSpotData[i][2],moRightPsycho[2])]

        for ii in controlSpots:
            if spotX == spots.iloc[ii][0] and spotY == spots.iloc[ii][1]:
                controlSpotPsycho[0] = psychoSpotData[ii][0]
                controlSpotPsycho[1] = [i+j for i,j in zip(psychoSpotData[ii][1],controlSpotPsycho[1])]
                controlSpotPsycho[2] = [np.nanmean([i,j]) for i,j in zip(psychoSpotData[ii][2],controlSpotPsycho[2])] 
                       
        


        ## plotting psychometrics for different cortical groups
    lines = []
    print('Fitting psychometrics, {}/{} Done'.format(subIdx,len(data)))
    fig = plt.figure()
    dotColors = ['.b','.r','.g']
    lineColors = ['-b','-r','-g']
    plotCount = 0
    for psychoPlot in [visLeftPsycho, visRightPsycho, controlSpotPsycho]:
        fitParams = []
        fitLikes = []
        for repeat in range(10):
            params, L = psychofit.mle_fit_psycho(psychoPlot,
                                        P_model='erf_psycho_2gammas',
                                        parstart=np.array([0,20,.05,.05]),
                                        parmin=np.array([-5, 0., 0., 0.]),
                                        parmax=np.array([5, 100., 1, 1]),
                                        nfits=50)
            fitParams.append(params)
            fitLikes.append(L)
        # find the best params (with the lowest neg likelihood)
        params = fitParams[np.where(min(fitLikes))[0][0]]
        spotBias[i][subIdx] = params[0]
        spotSlope[i][subIdx] = params[1]
        spotLapseLow[i][subIdx] = params[2]
        spotLapseHigh[i][subIdx] = params[3]

        spotFits.append(params)
        #plot the psychometrics
        fitx = np.linspace(-1, 1, 100)
        fity = psychofit.erf_psycho_2gammas(params, fitx)

        line = plt.plot(fitx, fity, lineColors[plotCount])
        lines.append(line)
        plt.plot(psychoPlot[0], np.array(psychoPlot[2]), dotColors[plotCount])
        plotCount+=1
    ## Formatting the psychometric figure    
    LvisLine = mpl.lines.Line2D([],[],color='blue',marker='.',label='L Vis Off')
    RvisLine = mpl.lines.Line2D([],[],color='red',marker='.',label='R Vis Off')
    cLine = mpl.lines.Line2D([],[],color='green',marker='.',label='Control Spots')
    plt.legend(handles=[LvisLine,RvisLine,cLine])
    if len(np.unique(subject['subject'])) == 1:
        fig.suptitle(np.unique(subject['subject'])[0])
    else:
        fig.suptitle('{} animals'.format(len(np.unique(subject['subject']))))


    lines = []
    print('Fitting psychometrics, {}/{} Done'.format(subIdx,len(data)))
    fig = plt.figure()
    dotColors = ['.b','.r','.g']
    lineColors = ['-b','-r','-g']
    plotCount = 0
    for psychoPlot in [moLeftPsycho, moRightPsycho, controlSpotPsycho]:
        fitParams = []
        fitLikes = []
        for repeat in range(10):
            params, L = psychofit.mle_fit_psycho(psychoPlot,
                                        P_model='erf_psycho_2gammas',
                                        parstart=np.array([0,20,.05,.05]),
                                        parmin=np.array([-5, 0., 0., 0.]),
                                        parmax=np.array([5, 100., 1, 1]),
                                        nfits=50)
            fitParams.append(params)
            fitLikes.append(L)
        # find the best params (with the lowest neg likelihood)
        params = fitParams[np.where(min(fitLikes))[0][0]]
        spotBias[i][subIdx] = params[0]
        spotSlope[i][subIdx] = params[1]
        spotLapseLow[i][subIdx] = params[2]
        spotLapseHigh[i][subIdx] = params[3]

        spotFits.append(params)
        #plot the psychometrics
        fitx = np.linspace(-1, 1, 100)
        fity = psychofit.erf_psycho_2gammas(params, fitx)

        line = plt.plot(fitx, fity, lineColors[plotCount])
        lines.append(line)
        plt.plot(psychoPlot[0], np.array(psychoPlot[2]), dotColors[plotCount])
        plotCount+=1
    ## Formatting the psychometric figure    
    LvisLine = mpl.lines.Line2D([],[],color='blue',marker='.',label='L Mo Off')
    RvisLine = mpl.lines.Line2D([],[],color='red',marker='.',label='R Mo Off')
    cLine = mpl.lines.Line2D([],[],color='green',marker='.',label='Control Spots')
    plt.legend(handles=[LvisLine,RvisLine,cLine])
    if len(np.unique(subject['subject'])) == 1:
        fig.suptitle(np.unique(subject['subject'])[0])
    else:
        fig.suptitle('{} animals'.format(len(np.unique(subject['subject']))))


    pSizes = []
    useToPlot = goRight
    plotLabels = ['Percent CCW', 'Percent CW', 'Percent No Go', 'Response Time', 'Percent Correct']
    for p in pVals:
        # Bonferroni correction, dev by num spots
        if p[useToPlot] <= .0001 / len(spots):
            pSizes.append(300)
        elif p[useToPlot] <= .001 / len(spots):
            pSizes.append(200)
        elif p[useToPlot] <= .01 / len(spots):
            pSizes.append(100)
        else:
            pSizes.append(5)

    fig = plt.figure()
    if len(np.unique(subject['subject'])) == 1:
        fig.suptitle(np.unique(subject['subject'])[0])
    else:
        fig.suptitle('{} animals'.format(len(np.unique(subject['subject']))))
    fig = plt.subplot2grid((5, 5), (0, 0), colspan=4, rowspan=4)
    plt.imshow(allenOutline, cmap="gray")
    plt.plot(bregma[0],bregma[1],'xk')
    allenSpotsX = (spots.iloc[:, 0] * 1 / pixelsize) + bregma[0]
    allenSpotsY = ((spots.iloc[:, 1]) * -1 / pixelsize) + bregma[1]

    useHue = np.array([mean[useToPlot] for mean in means])*100
    #cmap = sns.color_palette("RdBu_r", len(np.unique(useHue)))
    cmap = mpl.cm.seismic
    norm = mpl.colors.Normalize(vmin=0, vmax=100)
    textColor = 'w'
    if useToPlot == RT:
        useHue = np.array([mean[useToPlot] for mean in means])
        #cmap = sns.cubehelix_palette(len(np.unique(useHue)), start=1, rot=0, dark=0, light=.95)
        cmap = mpl.cm.Blues
        norm = mpl.colors.Normalize(vmin=0, vmax=0.5)
        textColor = 'k'
    elif useToPlot == correct:
        #cmap = sns.cubehelix_palette(len(np.unique(useHue)), start=1, rot=0, dark=0, light=.95, reverse=True)
        cmap = mpl.cm.Reds
        textColor = 'k'
    elif useToPlot == noGo:
        cmap = mpl.cm.Reds
        textColor = 'k'

    ax = sns.scatterplot(x=allenSpotsX, y=allenSpotsY, size=pSizes, sizes=(20, 400),
                         hue=useHue, palette=cmap, legend=False, edgecolor='k', hue_norm = norm)

    for x, y, p, h in zip(allenSpotsX,allenSpotsY,pSizes,useHue):
        if p == 300:
            plt.text(x,y,str(int(h)),horizontalalignment='center', verticalalignment='center', color=textColor)

    plt.text(max(allenSpotsX) + 80, 180, 'n = {} trials'.format(numTrials))
    if len(plotContrasts) == 5:
        plt.text(max(allenSpotsX) + 80, 200, 'All Contrasts')
    else:
        plt.text(max(allenSpotsX) + 80, 200, 'Contrast set: ' + str(plotContrasts))
    ax = plt.gca()
    ax.set_facecolor('w')
    plt.axis('off')
    maxColor = round(max(useHue), 2)
    minColor = round(min(useHue), 2)
    ax1 = plt.subplot2grid((5, 5), (4, 0), colspan=4, rowspan=1)
    colorBar = pd.DataFrame(np.sort(useHue[np.newaxis, :]), index=[plotLabels[useToPlot]],
                            columns=[str(int(round(hue * 100, 0))) for hue in np.sort(useHue)])
    #sns.heatmap(colorBar, cmap=cmap, cbar=False, ax=ax1, xticklabels=11, yticklabels=False)
    cb = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, orientation='horizontal', norm=norm)
    cb.set_label(plotLabels[useToPlot])
    
    legend_el = [mpl.lines.Line2D([0], [0], marker='o', color='w', label='p < .0001', markerfacecolor='k',
                        markersize=23),
                 mpl.lines.Line2D([0], [0], marker='o', color='w', label='p < .001', markerfacecolor='k',
                        markersize=18),
                 mpl.lines.Line2D([0], [0], marker='o', color='w', label='p < .01', markerfacecolor='k',
                        markersize=13),
                 mpl.lines.Line2D([0], [0], marker='o', color='w', label='p > .01', markerfacecolor='k',
                        markersize=6)]

    fig.legend(handles=legend_el,  loc='lower right', bbox_to_anchor=(1.5,.25), frameon=False,
              facecolor='w', columnspacing=10)
    subIdx+=1




def plot_psychos_from_spots(spotlists, psychoSpotData, spots, labels):
    """
        Function to plot three psychometrics on a single plot, with the dots that they were 
        generated with. usually I'll use this with a group of left spots, a group of right spots
        and a group of control spots.

        Inputs:
        spotlists: a list of length 3 that has the spots 

    """
    Psycho1 = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]
    Psycho2 = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]
    controlSpotPsycho = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]

    for i in range(len(spots)):
        spotX = spots.iloc[i, [0]][0]
        spotY = spots.iloc[i, [1]][0]
        for leftSpots in range(len(visLeftSpots)):
            if spotX == visLeftSpots[leftSpots][0] and spotY == visLeftSpots[leftSpots][1]:
                Psycho1[0] = psychoSpotData[i][0]
                Psycho1[1] = [i+j for i,j in zip(psychoSpotData[i][1],Psycho1[1])]
                Psycho1[2] = [np.nanmean([i,j]) for i,j in zip(psychoSpotData[i][2],Psycho1[2])]
        for rightSpots in range(len(visRightSpots)):
            if spotX == visRightSpots[rightSpots][0] and spotY == visRightSpots[rightSpots][1]:
                Psycho2[0] = psychoSpotData[i][0]
                Psycho2[1] = [i+j for i,j in zip(psychoSpotData[i][1],Psycho2[1])]
                Psycho2[2] = [np.nanmean([i,j]) for i,j in zip(psychoSpotData[i][2],Psycho2[2])]
        for ii in controlSpots:
            if spotX == spots.iloc[ii][0] and spotY == spots.iloc[ii][1]:
                controlSpotPsycho[0] = psychoSpotData[ii][0]
                controlSpotPsycho[1] = [i+j for i,j in zip(psychoSpotData[ii][1],controlSpotPsycho[1])]
                controlSpotPsycho[2] = [np.nanmean([i,j]) for i,j in zip(psychoSpotData[ii][2],controlSpotPsycho[2])] 
                       
        


        ## plotting psychometrics for different cortical groups
    lines = []
    print('Fitting psychometrics, {}/{} Done'.format(subIdx,len(data)))
    fig = plt.figure()
    dotColors = ['.b','.r','.g']
    lineColors = ['-b','-r','-g']
    plotCount = 0
    for psychoPlot in [Psycho1, Psycho2, controlSpotPsycho]:
        fitParams = []
        fitLikes = []
        for repeat in range(10):
            params, L = psychofit.mle_fit_psycho(psychoPlot,
                                        P_model='erf_psycho_2gammas',
                                        parstart=np.array([0,20,.05,.05]),
                                        parmin=np.array([-5, 0., 0., 0.]),
                                        parmax=np.array([5, 100., 1, 1]),
                                        nfits=50)
            fitParams.append(params)
            fitLikes.append(L)
        # find the best params (with the lowest neg likelihood)
        params = fitParams[np.where(min(fitLikes))[0][0]]
        spotBias[i][subIdx] = params[0]
        spotSlope[i][subIdx] = params[1]
        spotLapseLow[i][subIdx] = params[2]
        spotLapseHigh[i][subIdx] = params[3]

        spotFits.append(params)
        #plot the psychometrics
        fitx = np.linspace(-1, 1, 100)
        fity = psychofit.erf_psycho_2gammas(params, fitx)

        line = plt.plot(fitx, fity, lineColors[plotCount])
        lines.append(line)
        plt.plot(psychoPlot[0], np.array(psychoPlot[2]), dotColors[plotCount])
        plotCount+=1
    ## Formatting the psychometric figure    
    LvisLine = mpl.lines.Line2D([],[],color='blue',marker='.',label=labels[0])
    RvisLine = mpl.lines.Line2D([],[],color='red',marker='.',label=labels[1])
    cLine = mpl.lines.Line2D([],[],color='green',marker='.',label=labels[2])
    plt.legend(handles=[LvisLine,RvisLine,cLine])