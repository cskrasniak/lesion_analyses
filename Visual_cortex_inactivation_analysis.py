# Visual cortex inactivation analysis, to analyze data from biasedVisOffChoiceWorld, inactivating
# visual cortex to see if mice can perform without activity there

import numpy as np
import pandas as pd
from ibl_pipeline.utils import psychofit
import query_scan_mice as query
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil
from scipy import stats
from matplotlib.lines import Line2D
import datetime


bregma = [228.5, 190]
pixelsize = .025  # mm
allenOutline = np.load(r'F:\allen_dorsal_outline')
all_data = query.align_laser2behavior(['CSK-scan-005'])
dates = np.unique(all_data[0]['session_start_date'])
data = all_data[0]
data = data.reset_index()
data.rename(columns = {'laserPosX':'laserOn'}, inplace = True)
laserProb = data['laserOn'].rolling(50).mean()
#laserProb[0:24] = [.3]*24

data['trial_feedback_type'][data['trial_feedback_type'] == -1] = 0
percentCorrect = data['trial_feedback_type'].rolling(50).mean()
lateData = data[data['session_start_date'] == max(data['session_start_date'])]

laserOn = data['laserOn'] == 1 # find when laser is on
afterLaser = np.concatenate((np.array([False]),laserOn)) # add a laser Off to the beginning, this gives an array with if the laser was on or off for the last trial
afterLaser = afterLaser[:-1] # remove the extra element at the end
data['last_trial_LaserOn'] = afterLaser
data['afterLaserOff'] = afterLaser == False
np.mean(data['trial_feedback_type'][afterLaser])
np.mean(data['trial_feedback_type'][afterLaser == False]) # trials where the laser was not on after the last trial

plt.figure()
sns.barplot(y='trial_feedback_type', data=data, x='last_trial_LaserOn', hue='session_start_date')


plt.figure()
plt.plot(data.index,laserProb)
plt.plot(data.index,percentCorrect)
plt.plot(data.index,np.ones((len(data.index),1))*.5,'k' )
plt.legend(['Laser Probability','Correct Probability'])
plt.ylabel('Probability')
plt.xlabel('Trials')

lowerConData = data[abs(data['signed_contrast'] < 1)]
plt.figure()
sns.barplot(y='trial_feedback_type', data=lowerConData, x='laserOn', hue='session_start_date')
