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
import platform

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
        2019-12 CK added cumulative normal function
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
    ii = np.isfinite(data[2,:])

    likelihoods = np.zeros(nfits,)
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

    dispatcher={'psych_cumNorm': psych_cumNorm, 'erf_psycho_2gammas' : erf_psycho_2gammas}
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
def erf_psycho_2gammas(pars, xx):
    """
    erf function from 0 to 1, with two lapse rates.

    Args:
        pars: Model parameters [threshold, slope, gamma].
        xx: vector of stim levels (%)

    Returns:
        ff: A vector of length xx
        
    Examples:
        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> xx = np.arange(-50,50)
        >>> ff = erf_psycho_2gammas(np.array(-10., 10., 0.2, 0.),xx)
        >>> plt.plot(xx,ff)
        
    Raises: 
        ValueError: pars must be a vector of length 3
        ValueError: each of the three parameters must be scalar
        TypeError: pars must be a list or numpy array
        
    Information:
        2000    MC wrote it
        2018-08 MW ported to Python
    """
    # Validate input
    if isinstance(pars, (list, tuple)):
        pars = np.array(pars)
    elif not isinstance(pars, np.ndarray):
        raise TypeError('pars must be a list or numpy array')

    if pars.shape[0] != 4:
        raise ValueError('pars must be a vector of length 4')
    threshold	= pars[0]
    slope	= pars[1]
    gamma1	= pars[2]
    gamma2	= pars[3]

    if (threshold.size!=1) or (slope.size!=1) or (gamma1.size!=1) or (gamma2.size!=1):
        ValueError('each of the three parameters must be scalar')
    
    return gamma1 + (1 - gamma1 - gamma2) * (erf( (xx-threshold)/slope ) + 1 )/2
def psych_cumNorm(params, x):
    '''
    Psychometric function based on the cumulative normal distribution. inputs and returns a param
    vector [mu, invSigma, gamma, lambda] or [category boundary, inverse slope, low lapse, high lapse]
    '''
    mu = params[0]  # mean of the gaussian, x value for chance level
    invSigma = params[1]  # inverse variance of gaussian, aka inverse slope of psycho
    gamma = params[2]  # lapse Low
    lamb = params[3]  # lapse High
    y = gamma + (1 - gamma - lamb) * .5 * erfc(-invSigma *(x-mu) / np.sqrt(2))
    return y

def plot_from_spots(spotlist, psychoSpotData, spots, color, ax, plotType='psycho'):
    """
        Function to plot three psychometrics on a single plot, with the dots that they were 
        generated with. usually I'll use this with a group of left spots, a group of right spots
        and a group of control spots.

        Inputs:
        spotlists: an array containing the spots that you want to group together to plot the 
        psychometric for. a n by 2 numpy array that gives [[x1,y1],[x2,y2],[xn,yn]]
        psychoSpotData: a list of len(numSpots) each with a list lenght 4 containing 1, an array of
        signed contrasts, 2, a list of number of presentations for that contrast, 3, a list of 
        arrays containing a bool for if choice was CCW (list len num contrast, array len num
        presentations) 4, a list of the same form as above with the reaction times
        spots: a dataframe of length num spots, column 0: laserPosX, column 1: laserPosY, column3, 
        count
        color: the color for the line you want to plot, string eg 'b'
        ax: the matplotlib axis to plot onto
        plotType: defualt is psycho for plotting psychometric curves, other option is 'chrono'

    """
    if plotType == 'psycho':
        psycho = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.array([])]*len(psychoSpotData[0][1])]
        controlSpotPsycho = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.nan for i in range(len(psychoSpotData[0][1]))]]
        tempPsych = [np.array([])]*len(psychoSpotData[0][1])
        for i in range(len(spots)):
            spotX = spots.iloc[i, [0]][0]
            spotY = spots.iloc[i, [1]][0]
            for spot in range(len(spotlist)):
                if spotX == spotlist[spot][0] and spotY == spotlist[spot][1]:
                    psycho[0] = psychoSpotData[i][0]
                    psycho[1] = [temp+j for temp,j in zip(psychoSpotData[i][1],psycho[1])]
                    for contrast in range(len(psycho[0])):
                        con = int(contrast)
                        tempPsych[con] = np.append(tempPsych[con], psychoSpotData[i][2][con])

        sems = []
        for c in range(len(tempPsych)):
            psycho[2][c] = np.nanmean(tempPsych[c])
            sems.append(stats.sem(tempPsych[c]))

        ## Bootstrap confidence intervals
        nboots = 10
        bootFits = pd.DataFrame(columns=['threshold','slope','gamma','lambda'], index=range(nboots))
        bootData = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.array([])]*len(psychoSpotData[0][1])]
        bootData[0] = psycho[0]
        cnt=0
        print('bootstrapping errorbars...', sep=' ', end='')
        for i in range(nboots):
            if not(cnt % 5):
                print(int(cnt/nboots*100), sep=' ', end='%,', flush=True)
            for j in range(len(tempPsych)):
                bootData[2][j] = np.random.choice(tempPsych[j],size=int(len(tempPsych[j])/1.25), replace=True)                
                bootData[1][j] = len(bootData[2][j])
                bootData[2][j] = np.mean(bootData[2][j])
            fitParams = []
            fitLikes = []
            for repeat in range(5):
                parStart = np.array([-5+np.random.rand()*10,0+np.random.rand()*100, 0+np.random.rand(), 0+np.random.rand()])
                pars, L = mle_fit_psycho(bootData, P_model='erf_psycho_2gammas',
                                                            parstart=np.array([0, 50, .5, .5]),
                                                            parmin=np.array([-5, 0., 0., 0.]),
                                                            parmax=np.array([5, 100., 1, 1]),
                                                            nfits=2)
                
                fitParams.append(pars)
                fitLikes.append(L)
            cnt+=1
            bootFits.iloc[i] = fitParams[np.where(min(fitLikes))[0][0]]
            

        a = .05
        CIs = []
        for i in bootFits.columns:
            CIs.append([np.percentile(bootFits[i],100-a/2), np.percentile(bootFits[i],a/2,)])
            


            ## plotting psychometrics for different cortical groups
        lines = []
        
        fitParams = []
        fitLikes = []
        for repeat in range(10):
            parStart = np.array([-5+np.random.rand()*10,0+np.random.rand()*100, 0+np.random.rand(), 0+np.random.rand()])
            params, L = psychofit.mle_fit_psycho(psycho,
                                        P_model='erf_psycho_2gammas',
                                        parstart=parStart,
                                        parmin=np.array([-5, 0., 0., 0.]),
                                        parmax=np.array([5, 100., 1, 1]),
                                        nfits=25)
            fitParams.append(params)
            fitLikes.append(L)
            # find the best params (with the lowest neg likelihood)
        params = fitParams[np.where(min(fitLikes))[0][0]]

        spotFits.append(params)
        #plot the psychometrics
        fitx = np.linspace(-1, 1, 100)
        fity = psychofit.erf_psycho_2gammas(params, fitx)

        
        line = ax.plot(fitx, fity, color=color)
        lines.append(line)
        ax.errorbar(psycho[0], np.array(psycho[2]),yerr=sems, color=color, marker='.',ms=4, ls='')
        plt.ylim(-0.1,1.1)
        plt.xlim(-1.1,1.1)
    
    elif plotType == 'chrono':
        chrono = [[],[0 for i in range(len(psychoSpotData[0][1]))],[np.array([])]*len(psychoSpotData[0][1])]
        tempChrono = [np.array([])]*len(psychoSpotData[0][1])
        for i in range(len(spots)):
            spotX = spots.iloc[i, [0]][0]
            spotY = spots.iloc[i, [1]][0]
            
            for spot in range(len(spotlist)):
                if spotX == spotlist[spot][0] and spotY == spotlist[spot][1]:
                    chrono[0] = psychoSpotData[i][0]
                    chrono[1] = [temp+j for temp,j in zip(psychoSpotData[i][1],chrono[1])]
                    for contrast in range(len(chrono[0])):
                        con = int(contrast)
                        tempChrono[con] = np.append(tempChrono[con], psychoSpotData[i][3][con])
        sems = []
        for c in range(len(tempChrono)):
            chrono[2][c] = np.nanmedian(tempChrono[c])
            sems.append(stats.sem(tempChrono[c]))
        ax.errorbar(chrono[0], np.array(chrono[2]), yerr=sems, color=color, marker='.',ms=4)
        ax.set_ylim(.15,.5)
# the params I use here are RT at 0 contrast, 'RT bias' which is pairwise RT left- RT right, and 
# peakiness which is the ratio of max RT to the average of the two 100% RTs
        p2 = (chrono[2][0] - chrono[2][-1]) + (chrono[2][1] - chrono[2][-2] + (chrono[2][2] - chrono[2][-3]) + (chrono[2][3] - chrono[2][-4]))
        params = [chrono[2][4], p2, max(chrono[2])/np.mean([chrono[2][0],chrono[2][-1]])] 
        CIs = None
    else:
        raise Exception("This is not a supported plot type, choose 'psycho' or 'chrono'")
    return params, CIs

def err_bar_from_CIs(allCIs, x, y, ax):

    for i in range(len(allCIs)):

        ax.axvline(x=allCIs[x][i], ymin=allCIs[y][i][1], ymax=allCIs[y][i][0],lw=1,color='k')
        

################################# Start of script #########################################

bregma = [228.5, 190]
pixelsize = .025  # mm
if platform.system() == 'Darwin':
    allenOutline = np.load('/Users/ckrasnia/Desktop/Zador_Lab/scanData/allen_dorsal_outline')
else:
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
subList = []

plotContrasts = [1, .25, .125, .0625, 0]  # 1, .25, .125, list of contrasts to use in spot plot

visLeftSpots = np.array([[-2.5,-1.5],[-3.5,-1.5],[-2.5,-2.5],[-3.5,-2.5]])
visRightSpots = np.array([[2.5,-1.5],[3.5,-1.5],[2.5,-2.5],[3.5,-2.5]])
visPsychoParams = [[],[],[],[],[]]
visPsychoCIs = [[],[],[],[],[]]
visChronoParams = [[],[],[],[],[]]

moLeftSpots = np.array([[-1.5,.5],[-1.5,1.5],[-1.5,2.5],[-2.5,1.5]])
moRightSpots = np.array([[1.5,.5],[1.5,1.5],[1.5,2.5],[2.5,1.5]])
moPsychoParams = [[],[],[],[],[]]
moPsychoCIs = [[],[],[],[],[]]
moChronoParams = [[],[],[],[],[]]

for subject in data:

    spots = subject.groupby(['laserPosX', 'laserPosY']).size().reset_index().rename(
        columns={0: 'count'})
    controlSpots = list(spots[spots['laserPosY'] <-6].index)  # the spots that are very caudal
    contrastSet = np.unique(subject['signed_contrast'])
    subject['CCW'] = subject['trial_response_choice'] == 'CCW'
    subject = subject[subject.trial_response_choice != 'No Go']  # removes no go responses, toggle on and off
    subject.reset_index(drop=True)
    psychoSpotData = [[contrastSet, [], [], []] for spot in range(len(spots))]
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
    
    for i in range(len(spots)):
        spot = spots.iloc[i, [0, 1]]
        tempX = subject[subject['laserPosX'] == spot[0]]
        tempSpot = tempX[tempX['laserPosY'] == spot[1]]

        for contrast in contrastSet:
            byContrast = tempSpot[tempSpot['signed_contrast'] == contrast]
            psychoSpotData[i][1].append(len(byContrast))
            # psychoSpotData[i][2].append(np.mean(byContrast['trial_feedback_type']))
            psychoSpotData[i][2].append(np.array(byContrast['CCW']))
            psychoSpotData[i][3].append(np.array(byContrast['trial_response_time'] - byContrast[
                    'trial_go_cue_trigger_time']))
            if abs(contrast) in plotContrasts:
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


        ## plotting psychometrics for different cortical groups
    # plot formatting
    print('Fitting psychometrics')
    fig,axs = plt.subplots(nrows=2,ncols=2)
    ax1 = axs[0,0]
    ax1.set_title('Visual Cortex')
    ax1.set_ylabel('Fraction Choose Left')
    ax2 = axs[1,0]
    ax2.set_ylabel('Median RT (s)')
    ax2.set_xlabel('Signed Contrast')
    colors = ['r','b','g']
    # psychos and chronos for visual cortex
    plotCount = 0
    controlCoords = np.array([[spots.iloc[i,0],spots.iloc[i,1]] for i in controlSpots])
    
    for psychoSpots in [visLeftSpots, visRightSpots, controlCoords]:
        vPP, vPC = plot_from_spots(psychoSpots, psychoSpotData, spots, colors[plotCount], ax1, plotType='psycho')
        visPsychoParams[subIdx].append(vPP)
        visPsychoCIs[subIdx].append(vPC)
        vCP, vCC = plot_from_spots(psychoSpots, psychoSpotData, spots, colors[plotCount], ax2, plotType='chrono')
        visChronoParams[subIdx].append(vCP)
        plotCount+=1
    # psychos and chronos for motor cortex
    ax1 = axs[0,1]
    ax1.set_title('Motor Cortex')
    ax2 = axs[1,1]
    ax2.set_xlabel('Signed Contrast')
    plotCount = 0
    for psychoSpots in [moLeftSpots, moRightSpots, controlCoords]:
        mPP, mPC = plot_from_spots(psychoSpots, psychoSpotData, spots, colors[plotCount], ax1, plotType='psycho')
        moPsychoParams[subIdx].append(mPP)
        moPsychoCIs[subIdx].append(mPC)
        mCP, mCC = plot_from_spots(psychoSpots, psychoSpotData, spots, colors[plotCount], ax2, plotType='chrono')
        moChronoParams[subIdx].append(mCP)
        plotCount+=1
    
    ## making the figure legend and title  
    LvisLine = mpl.lines.Line2D([],[],color='red',marker='.',label='Laser Left')
    RvisLine = mpl.lines.Line2D([],[],color='blue',marker='.',label='Laser Right')
    cLine = mpl.lines.Line2D([],[],color='green',marker='.',label='Control')
    ax2.legend(handles=[LvisLine,RvisLine,cLine],loc='upper right')
    if len(np.unique(subject['subject'])) == 1:
        fig.suptitle(np.unique(subject['subject'])[0])
    else:
        fig.suptitle('{} animals'.format(len(np.unique(subject['subject']))))

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
    ## sorting Pvals into different categories
    pSizes = []
    useToPlot = goLeft
    plotLabels = ['Percent Leftward Choice', 'Percent Choose Right', 'Percent No Go', 'Response Time', 'Percent Correct']
    for p in pVals:
        # Bonferroni correction, dev by num spots
        if p[useToPlot] <= .0001 / len(spots):
            pSizes.append(300)
        elif p[useToPlot] <= .001 / len(spots):
            pSizes.append(175)
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
    plt.imshow(allenOutline, cmap="gray", interpolation='nearest')
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
        if p >= 210:
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
                        markersize=21),
                 mpl.lines.Line2D([0], [0], marker='o', color='w', label='p < .001', markerfacecolor='k',
                        markersize=17),
                 mpl.lines.Line2D([0], [0], marker='o', color='w', label='p < .01', markerfacecolor='k',
                        markersize=13),
                 mpl.lines.Line2D([0], [0], marker='o', color='w', label='p > .01', markerfacecolor='k',
                        markersize=6)]

    fig.legend(handles=legend_el,  loc='lower right', bbox_to_anchor=(1.5,.25), frameon=False,
              facecolor='w', columnspacing=10)
    subIdx+=1
    if len(np.unique(subject['subject'])) == 1:
        subList.append(np.unique(subject['subject'])[0])
    else:
        subList.append('{} animals'.format(len(np.unique(subject['subject']))))
plt.show(block=False)

## subject comparison analyses
allParams = pd.DataFrame(index=range(len(subList)*3), columns=['thresh', 'slope', 'gamma', 'lambda','subject', 'laserLocation'])
allCIs = pd.DataFrame(index=range(len(subList)*3), columns=['thresh', 'slope', 'gamma', 'lambda','subject', 'laserLocation'])
visCP = pd.DataFrame(index=range(len(subList)*3), columns=['RT0','RT bias', 'RT peakiness','subject', 'laserLocation'])
moCP = pd.DataFrame(index=range(len(subList)*3), columns=['RT0','RT bias', 'RT peakiness','subject', 'laserLocation'])
laserLocs = ['Left','Right','Control']
cnt = 0
for i in range(3):
    pars = visPsychoParams[:][i]
    for j in range(len(subList)):
        allParams.iloc[cnt,:4] = (visPsychoParams[j][i])
        allParams.iloc[cnt,4] = subList[j]
        allParams.iloc[cnt,5] = laserLocs[i]
        visCP.iloc[cnt,:3] = visChronoParams[j][i]
        visCP.iloc[cnt,3] = subList[j]
        visCP.iloc[cnt,4] = laserLocs[i]
        moCP.iloc[cnt,:3] = moChronoParams[j][i]
        moCP.iloc[cnt,3] = subList[j]
        moCP.iloc[cnt,4] = laserLocs[i]
        allCIs.iloc[cnt,:4] = visPsychoCIs[j][i]
        allCIs.iloc[cnt,4] = subList[j]
        allCIs.iloc[cnt,5] = laserLocs[i]
        cnt+=1

allParams.slope = 1/allParams.slope
pal = sns.color_palette('colorblind') 
fig, axs = plt.subplots(nrows=2,ncols=2)
fig.suptitle('Visual Cortex Inactivation')
ax1=axs[0,0]
sns.barplot(x='laserLocation', y='thresh', data=allParams,ax=ax1,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='thresh', data=allParams, hue='subject',ax=ax1)
sns.pointplot(x='laserLocation', y='thresh', data=allParams,hue='subject', ax=ax1, palette=pal)
# err_bar_from_CIs(allCIs, 'laserLocation', 'thresh', ax1)
ax1.get_legend().remove()
ax1.set_xlabel('')
ax1.set_ylabel('detection threshold\n(contrast fraction)')

ax2 = axs[0,1]
sns.barplot(x='laserLocation', y='slope', data=allParams,ax=ax2,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='slope', data=allParams, hue='subject',ax=ax2)
sns.pointplot(x='laserLocation', y='slope', data=allParams,hue='subject', ax=ax2, palette=pal)
# err_bar_from_CIs(allCIs, 'laserLocation', 'slope', ax2)
ax2.get_legend().remove()
ax2.set_xlabel('')
ax2.set_ylabel('slope\n(choose left/contrast)')

ax3 = axs[1,0]
sns.barplot(x='laserLocation', y='gamma', data=allParams,ax=ax3,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='gamma', data=allParams, hue='subject',ax=ax3)
# err_bar_from_CIs(allCIs, 'laserLocation', 'gamma', ax3)
sns.pointplot(x='laserLocation', y='gamma', data=allParams,hue='subject', ax=ax3, palette=pal)
ax3.get_legend().remove()
ax3.set_ylabel('right lapse\n(fraction contrast)')

ax4 = axs[1,1]
sns.barplot(x='laserLocation', y='lambda', data=allParams,ax=ax4,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='lambda', data=allParams, hue='subject',ax=ax4)
# err_bar_from_CIs(allCIs, 'laserLocation', 'lambda', ax4)
sns.pointplot(x='laserLocation', y='lambda', data=allParams,hue='subject', ax=ax4, palette=pal)
ax4.get_legend().remove()
ax4.set_ylabel('left lapse\n(fraction contrast)')

plt.show(block=False)

moAllParams = pd.DataFrame(index=range(len(subList)*3), columns=['thresh', 'slope', 'gamma', 'lambda','subject', 'laserLocation'])
moAllCIs = pd.DataFrame(index=range(len(subList)*3), columns=['thresh', 'slope', 'gamma', 'lambda','subject', 'laserLocation'])
cnt=0
for i in range(3):
    pars = visPsychoParams[:][i]
    for j in range(len(subList)):
        moAllParams.iloc[cnt,:4] = (moPsychoParams[j][i])
        moAllParams.iloc[cnt,4] = subList[j]
        moAllParams.iloc[cnt,5] = laserLocs[i]
        moAllCIs.iloc[cnt,:4] = moPsychoCIs[j][i]
        moAllCIs.iloc[cnt,4] = subList[j]
        moAllCIs.iloc[cnt,5] = laserLocs[i]
        cnt+=1
moAllParams.slope = 1/moAllParams.slope
fig, axs = plt.subplots(nrows=2,ncols=2)
fig.suptitle('Motor Cortex Inactivation')
ax1=axs[0,0]
sns.barplot(x='laserLocation', y='thresh', data=moAllParams,ax=ax1,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='thresh', data=moAllParams, hue='subject',ax=ax1)
sns.pointplot(x='laserLocation', y='thresh', data=moAllParams,hue='subject', ax=ax1, palette=pal)
# err_bar_from_CIs(allCIs, 'laserLocation', 'thresh', ax1)
ax1.get_legend().remove()
ax1.set_xlabel('')
ax1.set_ylabel('detection threshold\n(contrast fraction)')

ax2 = axs[0,1]
sns.barplot(x='laserLocation', y='slope', data=moAllParams,ax=ax2,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='slope', data=moAllParams, hue='subject',ax=ax2)
sns.pointplot(x='laserLocation', y='slope', data=moAllParams,hue='subject', ax=ax2, palette=pal)
# err_bar_from_CIs(allCIs, 'laserLocation', 'slope', ax2)
ax2.get_legend().remove()
ax2.set_xlabel('')
ax2.set_ylabel('slope\n(choose left/contrast)')

ax3 = axs[1,0]
sns.barplot(x='laserLocation', y='gamma', data=moAllParams,ax=ax3,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='gamma', data=moAllParams, hue='subject',ax=ax3)
# err_bar_from_CIs(allCIs, 'laserLocation', 'gamma', ax3)
sns.pointplot(x='laserLocation', y='gamma', data=moAllParams,hue='subject', ax=ax3, palette=pal)
ax3.get_legend().remove()
ax3.set_ylabel('right lapse\n(fraction contrast)')

ax4 = axs[1,1]
sns.barplot(x='laserLocation', y='lambda', data=moAllParams,ax=ax4,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='lambda', data=moAllParams, hue='subject',ax=ax4)
# err_bar_from_CIs(allCIs, 'laserLocation', 'lambda', ax4)
sns.pointplot(x='laserLocation', y='lambda', data=moAllParams,hue='subject', ax=ax4, palette=pal)
ax4.get_legend().remove()
ax4.set_ylabel('left lapse\n(fraction contrast)')

plt.show(block=False)


fig, axs = plt.subplots(nrows=2,ncols=2)
fig.suptitle('Visual Cortex Inactivation')
ax1=axs[0,0]
sns.barplot(x='laserLocation', y='RT0', data=visCP,ax=ax1,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='RT0', data=visCP, hue='subject',ax=ax1)
sns.pointplot(x='laserLocation', y='RT0', data=visCP,hue='subject', ax=ax1, palette=pal)
# err_bar_from_CIs(allCIs, 'laserLocation', 'RT0', ax1)
ax1.get_legend().remove()
ax1.set_xlabel('')
ax1.set_ylabel('Reaction time at 0 contrast\n(s)')

ax2 = axs[1,0]
sns.barplot(x='laserLocation', y='RT bias', data=visCP,ax=ax2,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='RT bias', data=visCP, hue='subject',ax=ax2)
sns.pointplot(x='laserLocation', y='RT bias', data=visCP,hue='subject', ax=ax2, palette=pal)
# err_bar_from_CIs(allCIs, 'laserLocation', 'RT bias', ax2)
ax2.get_legend().remove()
ax2.set_xlabel('')
ax2.set_ylabel('Difference between left and right\n reaction times (s)')

ax4 = axs[1,1]
sns.barplot(x='laserLocation', y='RT peakiness', data=visCP,ax=ax4,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='RT peakiness', data=visCP, hue='subject',ax=ax4)
# err_bar_from_CIs(allCIs, 'laserLocation', 'RT peakiness', ax4)
sns.pointplot(x='laserLocation', y='RT peakiness', data=visCP,hue='subject', ax=ax4, palette=pal)
ax4.get_legend().remove()
ax4.set_ylabel('ratio of RT for high\nand low contrasts')

plt.show(block=False)


fig, axs = plt.subplots(nrows=2,ncols=2)
fig.suptitle('Motor Cortex Inactivation')
ax1=axs[0,0]
sns.barplot(x='laserLocation', y='RT0', data=moCP,ax=ax1,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='RT0', data=moCP, hue='subject',ax=ax1)
sns.pointplot(x='laserLocation', y='RT0', data=moCP,hue='subject', ax=ax1, palette=pal)
# err_bar_from_CIs(allCIs, 'laserLocation', 'RT0', ax1)
ax1.get_legend().remove()
ax1.set_xlabel('')
ax1.set_ylabel('Reaction time at 0 contrast\n(s)')

ax2 = axs[1,0]
sns.barplot(x='laserLocation', y='RT bias', data=moCP,ax=ax2,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='RT bias', data=moCP, hue='subject',ax=ax2)
sns.pointplot(x='laserLocation', y='RT bias', data=moCP,hue='subject', ax=ax2, palette=pal)
# err_bar_from_CIs(allCIs, 'laserLocation', 'RT bias', ax2)
ax2.get_legend().remove()
ax2.set_xlabel('')
ax2.set_ylabel('Difference between left and right\n reaction times (s)')

ax4 = axs[1,1]
sns.barplot(x='laserLocation', y='RT peakiness', data=moCP,ax=ax4,ci=None, palette='Greys')
# sns.scatterplot(x='laserLocation', y='RT peakiness', data=moCP, hue='subject',ax=ax4)
# err_bar_from_CIs(allCIs, 'laserLocation', 'RT peakiness', ax4)
sns.pointplot(x='laserLocation', y='RT peakiness', data=moCP,hue='subject', ax=ax4, palette=pal)
ax4.get_legend().remove()
ax4.set_ylabel('ratio of RT for high\nand low contrasts')

plt.show(block=False)



vLeft = allParams.iloc[:4]
vRight = allParams.iloc[5:9]
vControl = allParams.iloc[10:14]
mLeft = moAllParams.iloc[:4]
mRight = moAllParams.iloc[5:9]
mControl = moAllParams.iloc[10:14]
anovaP =  pd.DataFrame(index=['vis','mo'], columns=['thresh', 'slope', 'gamma', 'lambda']) 
for i in anovaP.columns:
    f1, p1 =  scipy.stats.f_oneway(np.array(vLeft[i]), np.array(vRight[i]), np.array(vControl[i])) 
    f2, p2 = scipy.stats.f_oneway(np.array(mLeft[i]), np.array(mRight[i]), np.array(mControl[i]))
    anovaP[i]['vis'] = p1
    anovaP[i]['mo'] = p2
print('visResults:\n', anovaP)

vLeft = visCP.iloc[:4]
vRight = visCP.iloc[5:9]
vControl = visCP.iloc[10:14]
mLeft = moCP.iloc[:4]
mRight = moCP.iloc[5:9]
mControl = moCP.iloc[10:14]
anovaP =  pd.DataFrame(index=['vis','mo'], columns=['RT0', 'RT bias', 'RT peakiness']) 
for i in anovaP.columns:
    f1, p1 =  scipy.stats.f_oneway(np.array(vLeft[i]), np.array(vRight[i]), np.array(vControl[i])) 
    f2, p2 = scipy.stats.f_oneway(np.array(mLeft[i]), np.array(mRight[i]), np.array(mControl[i]))
    anovaP[i]['vis'] = p1
    anovaP[i]['mo'] = p2
print('moResults:\n', anovaP)
