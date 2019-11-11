import numpy as np
import pandas as pd
from ibl_pipeline import behavior
from ibl_pipeline.utils import psychofit
from ibl_pipeline.analyses import behavior as behavior_analysis
import query_scan_mice as query
import seaborn as sns
import matplotlib.pyplot as plt
import os
from math import ceil
from scipy import stats
from matplotlib.lines import Line2D

bregma = [228.5, 190]
pixelsize = .025  #mm
allenOutline = np.load(r'F:\allen_dorsal_outline')
#data = query.align_laser2behavior(['CSK-scan-004', 'CSK-scan-005'])

for subject in bigData:
    shuffledData = subject.copy()
    shuffledData['laserPosX'] = np.random.permutation(shuffledData.loc[:,'laserPosX'])
    shuffledData['laserPosY'] = np.random.permutation(shuffledData.loc[:,'laserPosY'])

    spots = subject.groupby(['laserPosX', 'laserPosY']).size().reset_index().rename(
        columns={0: 'count'})
    contrastSet = np.unique(subject['signed_contrast'])
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

    f, axes = plt.subplots(13, 10, figsize=(14, 14), sharex=True, sharey=True)
    sns.despine(left=True, bottom=True)
    

    for i in range(len(spots)):
        spot = spots.iloc[i, [0, 1]]
        tempX = subject[subject['laserPosX'] == spot[0]]
        tempSpot = tempX[tempX['laserPosY'] == spot[1]]

        for contrast in contrastSet:
            byContrast = tempSpot[tempSpot['signed_contrast'] == contrast]
            psychoSpotData[i][1].append(len(byContrast))
            psychoSpotData[i][2].append(np.mean(byContrast['trial_feedback_type']))
            if abs(contrast) == .125 or abs(contrast) == .0625:
                spotData[i][goLeft].append(byContrast['trial_response_choice'] == -1)  # go left
                spotData[i][goRight].append(byContrast['trial_response_choice'] == 1)  # go right
                spotData[i][noGo].append(byContrast['trial_response_choice'] == 0)  # no go
                spotData[i][RT].append(byContrast['trial_response_time'] - byContrast['trial_go_cue_trigger_time'])
                byContrast[byContrast['trial_feedback_type'] == -1] = 0
                spotData[i][correct].append(byContrast['trial_feedback_type'])  # correct
        # Compiling lists
        spotData[i][goLeft] = [item for sublist in spotData[i][goLeft] for item in sublist]
        spotData[i][goRight] = [item for sublist in spotData[i][goRight] for item in sublist]
        spotData[i][noGo] = [item for sublist in spotData[i][noGo] for item in sublist]
        spotData[i][RT] = [item for sublist in spotData[i][RT] for item in sublist]
        spotData[i][correct] = [item for sublist in spotData[i][correct] for item in sublist]

        params, L = psychofit.mle_fit_psycho(psychoSpotData[i],
                                             P_model='erf_psycho_2gammas',
                                             parstart=np.array([.3, 5, .2, .2]),
                                             parmin=np.array([.005, 0., 0., 0.]),
                                             parmax=np.array([.5, 10., .45, .45]))
        spotFits.append(params)
        fitx = np.linspace(-100, 100, 100)
        fity = psychofit.erf_psycho_2gammas(params, fitx)
        sns.lineplot(fitx, fity, ax=axes[abs(int(spots.iloc[i, 1]) - 4),
                     int(ceil(spots.iloc[i, 0] + 4))], palette='gray')
        plt.axis('off')
        ax.tick_params(left=False, bottom=False)

    controlData[goLeft] = spotData[10][goLeft] + spotData[55][goLeft]
    controlData[goRight] = spotData[10][goRight] + spotData[55][goRight]
    controlData[noGo] = spotData[10][noGo] + spotData[55][noGo]
    controlData[RT] = spotData[10][RT] + spotData[55][RT]
    controlData[correct] = spotData[10][correct] + spotData[55][correct]

    for i in range(len(spots)):
        t, pVals[i][goLeft] = stats.ttest_ind(spotData[i][goLeft], controlData[goLeft])
        t, pVals[i][goRight] = stats.ttest_ind(spotData[i][goRight], controlData[goRight])
        t, pVals[i][noGo] = stats.ttest_ind(spotData[i][noGo], controlData[noGo])
        t, pVals[i][RT] = stats.ttest_ind(spotData[i][RT], controlData[RT])
        t, pVals[i][correct] = stats.ttest_ind(spotData[i][correct], controlData[correct])
        means[i][goLeft] = np.mean(spotData[i][goLeft])
        means[i][goRight] = np.mean(spotData[i][goRight])
        means[i][noGo] = np.mean(spotData[i][noGo])
        means[i][RT] = np.mean(spotData[i][RT])
        means[i][correct] = np.mean(spotData[i][correct])

    pSizes = []
    useToPlot = goRight
    plotLabels = [' Prob go left', ' Prob go right', ' Prob no go', ' Response Time', ' Correct']
    for p in pVals:
        if p[useToPlot] <= .0001:
            pSizes.append(400)
        elif p[useToPlot] <= .001:
            pSizes.append(250)
        elif p[useToPlot] <= .01:
            pSizes.append(50)
        else:
            pSizes.append(20)
    fig = plt.figure(figsize=(4.6, 7))
    plt.imshow(allenOutline, cmap="gray")
    allenSpotsX = (spots.iloc[:, 0] * 40) + bregma[0]
    allenSpotsY = ((spots.iloc[:, 1]) * -40) + bregma[1]

    useHue = np.array([mean[useToPlot] for mean in means]) - np.mean(controlData[useToPlot])
    cmap = sns.color_palette("RdBu_r", len(np.unique(useHue)))
    red = cmap[-1]
    blue = cmap[0]
    if useToPlot == RT:
        cmap = sns.cubehelix_palette(len(np.unique(useHue)), start=1, rot=0, dark=0, light=.95)
        red = cmap[-1]
        blue = cmap[0]
    elif useToPlot == correct:
        cmap = sns.cubehelix_palette(len(np.unique(useHue)), start=1, rot=0, dark=0, light=.95,reverse=True)
        red = cmap[-1]
        blue = cmap[0]
    ax = sns.scatterplot(x=allenSpotsX, y=allenSpotsY, size=pSizes, sizes=(20, 400),
                         hue=useHue, palette=cmap, legend=False, edgecolor='k')
    ax = plt.gca()
    ax.set_facecolor('w')
    plt.axis('off')
    maxColor = round(max([mean[useToPlot] for mean in means]) - np.mean(controlData[useToPlot]), 2)
    minColor = round(min([mean[useToPlot] for mean in means]) - np.mean(controlData[useToPlot]), 2)
    legend_el = [Line2D([0], [0], marker='o', color ='w', label='p < .0001', markerfacecolor='k', markersize = 20),
             Line2D([0], [0], marker='o', color='w', label='p < .001', markerfacecolor='k', markersize = 15),
             Line2D([0], [0], marker='o', color='w', label='p < .01', markerfacecolor='k', markersize = 8),
             Line2D([0], [0], marker='o', color='w', label='p > .01', markerfacecolor='k', markersize = 5),
             Line2D([0], [0], marker='o', color='w', label='+ ' + str(maxColor) + plotLabels[useToPlot], markerfacecolor = red, markersize = 15),
             Line2D([0], [0], marker='o', color='w', label=str(minColor) + plotLabels[useToPlot], markerfacecolor = blue, markersize = 15)]

    ax.legend(handles=legend_el, bbox_to_anchor=(0.5, 0.05), loc='center', frameon=False,
              facecolor='w', columnspacing=10)
