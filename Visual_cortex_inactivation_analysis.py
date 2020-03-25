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
all_data = query.align_laser2behavior(['CSK-scan-017'])
dates = np.unique(all_data[0]['session_start_date'])
data = all_data[0][all_data[0]['session_start_date'] >= dates[-2]]
data = data.reset_index()
data.rename(columns = {'laserPosX':'laserOn'}, inplace = True)
laserProb = data['laserOn'].rolling(50).mean()
#laserProb[0:24] = [.3]*24

data['trial_feedback_type'][data['trial_feedback_type'] == -1] = 0
percentCorrect = data['trial_feedback_type'].rolling(50).mean()

plt.figure()
plt.plot(data.index,laserProb)
plt.plot(data.index,percentCorrect)
plt.plot(data.index,np.ones((len(data.index),1))*.5,'k' )
plt.legend(['Laser Probability','Correct Probability'])
plt.ylabel('Probability')
plt.xlabel('Trials')