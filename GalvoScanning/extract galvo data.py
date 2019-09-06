# -*- coding: utf-8 -*-
"""
Created on Sat Aug 10 11:35:12 2019

@author: chris
"""

## Extract Data

import os
import pandas as pd
import numpy as np
from alf.extractors.biased_trials import extract_all
subject = "CSK-scan-002"
#def extractGalvoData(subject):
trainData = pd.DataFrame(columns= ['feedbackType','contrastLeft','contrastRight','probabilityLeft','choice','rewardVolume','feedback_times','stimOn_times','response_times','iti_dur','goCue_times','goCueTrigger_times','laserX','laserY','laserOn'])

dataPath = "F:\Subjects"

os.chdir(dataPath)    
days = os.listdir(subject)

for day in days:
    os.chdir(dataPath)  
    os.chdir(subject)
    runs = os.listdir(day)
    training = []
    laserData = []
    
    
    for run in runs:         
        os.chdir(os.path.join(dataPath,subject,day))
        try:
            currpath = os.getcwd()
            training.append(extract_all(run)) #get behavior data
        except:
            pass
        os.chdir(run)
        laserData.append(np.load("laserData")) #get laser data

    useMax = []
    for el in training:
        useMax.append(len(el['feedbackType']))
    training = training[useMax.index(max(useMax))] #only take the longest training session 
    laserData = laserData[useMax.index(max(useMax))] #same for laser
    #        if np.size(laserData,1)>3:
#            laserData[] = 'visLeft'
#            laserData[] = 'visRight'
#            laserData[] = 'offLeft'
#            laserData[] = 'offleft'
    if max(useMax) <= min(useMax)*10 and max(useMax) != min(useMax):
        print('WARNING, there are two files that are close to the same length in folder:')
        print(currpath)
    del training['session_path']
    del training['intervals']
    training = pd.DataFrame.from_dict(training)
    if len(training) < len(laserData):
        laserData = np.delete(laserData,-1,0)
    training['laserX']= laserData[:,0]
    training['laserY']=laserData[:,1]
    training['laserOn'] = laserData[:,2]
    trainData = trainData.append(training)
    #os.chdir(os.path.dirname(os.getcwd()))
