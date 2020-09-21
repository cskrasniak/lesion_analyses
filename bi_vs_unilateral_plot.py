# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 09:50:06 2020

@author: chris
"""

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

all_data = query.align_laser2behavior(['CSK-scan-019'])# data is from these two mice days 06/03/2020-06/05/2020
data = all_data[0]#.append(all_data[1])

data['trial_feedback_type'][data['trial_feedback_type'] == -1] = 0

# for i in range(len(data)):
#     if pd.isna(data.iloc[i]['laserLoc']):
#         data.at[i,'laserLoc'] = data.iloc[i]['laserOn']
        
droptrials = data[data['trial_response_choice'] == 'No Go'].index
data.drop(droptrials, inplace=True)
data['laserOn'] = np.nan
data['laserOn'][data['laserLoc'] == 'control'] = 0
data['laserOn'][data['laserLoc'] == 'left'] = 1
data['laserOn'][data['laserLoc'] == 'right'] = 2
data['laserOn'][data['laserLoc'] == 'bilateral'] = 3
data=data.sort_values(by='laserOn')

plt.figure()
sns.barplot(data=data, x='laserLoc', y='trial_feedback_type')
plt.xlabel('Laser Location')
plt.ylabel('Proportion Correct')
plt.ylim(0,1)