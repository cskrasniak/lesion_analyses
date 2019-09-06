# -*- coding: utf-8 -*-
"""
Created on Fri Aug  9 12:33:13 2019

@author: chris
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from os.path import join
import seaborn as sns
import datajoint as dj
from ibl_pipeline import subject, acquisition, action, behavior, reference
from ibl_pipeline.analyses import behavior as behavior_analysis
from alf.extractors.biased_trials import extract_all
from ibl_pipeline.utils import psychofit
from matplotlib.lines import Line2D
### Getting and aranging behavior data
#data = extract_all(r"C:\Users\chris\Downloads\2019-08-08\2019-08-08\001")
#del data['session_path']
#del data['intervals']
#data = pd.DataFrame.from_dict(data)


cnt=0
data = useData.copy()
subject = 'CSK-scan-002'
data = data.reset_index()
data['correct'] = data['feedbackType']
data['signedContrast'] = data['contrastLeft'].sub(data['contrastRight'], fill_value = 0)


  

for i in range(len(data)):
    cnt+=1
    
    if data['choice'][i] == 0 :
        data['correct'][i] = 0
        print(cnt)
    else:
        data['correct'][i] = data['feedbackType'][i]
    


## Formatting data to feed into psychofit
shuffledDF = data.copy()
shuffledDF['laserX'] = np.random.permutation(shuffledDF.loc[:,'laserX'])
shuffledDF['laserY'] = np.random.permutation(shuffledDF.loc[:,'laserY'])

#data = data.iloc[8000:]
#data = data.reset_index()
#shuffledDF = shuffledDF.iloc[8000:]
#shuffledDF = shuffledDF.reset_index()  

allPoints = np.load(r'F:\allPoints')
contrastSet = [-1., -.25, -.0625, 0, .0625, .25, 1]
signedContrasts = [ [] for i in range(len(allPoints)) ]
correct = [ [] for i in range(len(allPoints)) ]
response = [ [] for i in range(len(allPoints)) ]
signedContrasts_shuf = [ [] for i in range(len(allPoints)) ]
correct_shuf = [ [] for i in range(len(allPoints)) ]
response_shuf = [ [] for i in range(len(allPoints)) ]
psychFitters = [ np.array([[],[],[]]) for i in range(len(allPoints))]
psychFitters_all = np.array([contrastSet,np.zeros(len(contrastSet)),np.empty(len(contrastSet))])


#for k in range(len(contrastSet)):
#    psyCor = 0
#    psycnt = 0


for i in range(len(data)):    
    
    
#        if contrastSet[k] == data['signedContrast'].iloc[i]:
#            psychFitters_all[1,k]+=1
#            psycnt+=1
#            if data['feedbackType'].iloc[i] == 1:
#                psyCor +=1
#    psychFitters_all[2,k] = psyCor/psycnt
    
    
    for j in range(len(allPoints)):
        if data["laserOn"].iloc[i] == 1 and data["laserX"].iloc[i]==allPoints[j,0] and data['laserY'].iloc[i] == allPoints[j,1] and abs(data['signedContrast'][i])>=.25:
            correct[j].append(data['feedbackType'][i])
            response[j].append(data['choice'][i])
        if shuffledDF["laserOn"].iloc[i] == 1 and shuffledDF["laserX"].iloc[i]==allPoints[j,0] and shuffledDF['laserY'].iloc[i] == allPoints[j,1] and abs(shuffledDF['signedContrast'][i])>=.25 :
            signedContrasts_shuf[j].append(shuffledDF['signedContrast'][i])
            correct_shuf[j].append(shuffledDF['feedbackType'][i])
            response_shuf[j].append(shuffledDF['choice'][i])
        ## for fitting psychs
#        signedContrasts[j].append(data['signedContrast'][i])
#        for k in contrastSet:
#            if data['signedContrast'][i] == k:
#                psychFitters[j,1] += 1
      
#pars,l =  psychofit.mle_fit_psycho(psychFitters_all, parstart = np.array([5,5,5,5]),parmin = np.array([-1,0,0,0]),parmax =np.array([1,10,1,1]),  P_model = 'erf_psycho_2gammas')
#plotFit = psychofit.erf_psycho_2gammas(pars,np.linspace(-1,1,1000))
#plt.figure()
#sns.lineplot(x = np.linspace(-1,1,1000),y = plotFit)
        
        
p_vals = []
t_vals = []
responseBias = []
controls = [val for sublist in response[-2:] for val in sublist]
for i in range(len(response)):
    t,p = stats.ttest_ind(response[i],controls)
    p_vals.append(p)   
    t_vals.append(t)         
meanResponses = []
goRight = []
goLeft = []
noGo = []
for i in range(len(response)):
    meanResponses.append((np.mean(response[i])))
    goRight.append((response[i].count(-1)/len(response[i])) - controls.count(-1)/len(controls))
    goLeft.append((response[i].count(1)/len(response[i]))  - controls.count(-1)/len(controls))
    noGo.append((response[i].count(0)/len(response[i]))  - controls.count(0)/len(controls))
pSizes = []
for p in p_vals:
    if p <=.001:
       pSizes.append(400)
    elif p <= .01:
       pSizes.append(250)
    elif p <=.05:
       pSizes.append(50)
    else:
       pSizes.append(20)
plt.figure(figsize = (4.6,7))
sns.set_palette(sns.color_palette("RdBu_r"))   
red = sns.color_palette('RdBu_r',70)[-1]
blue = sns.color_palette('RdBu_r',70)[0]
ax = sns.scatterplot(x=allPoints[:,0], y = allPoints[:,1], size = pSizes, sizes =(20,400),  hue = noGo ,palette = "RdBu_r",legend= False,edgecolor = 'k')
ax = plt.gca()
ax.set_facecolor('w')
plt.axis('off')
plt.title('1 & 2, 0 Contrast')
legend_el = [Line2D([0],[0], marker = 'o', color = 'w',label = 'p < .001', markerfacecolor = 'k', markersize = 20),
             Line2D([0],[0], marker = 'o', color = 'w',label = 'p < .01', markerfacecolor = 'k', markersize = 15),
             Line2D([0],[0], marker = 'o', color = 'w',label = 'p < .05', markerfacecolor = 'k', markersize = 8),
             Line2D([0],[0], marker = 'o', color = 'w',label = 'p > .05', markerfacecolor = 'k', markersize = 5),
             Line2D([0],[0], marker = 'o', color = 'w', label = str(round(max(noGo),2)) +' \u0394 prob. no go', markerfacecolor = red, markersize = 15),
             Line2D([0],[0], marker = 'o', color = 'w', label = str(round(min(noGo),2)) +' \u0394 prob. no go', markerfacecolor = blue, markersize = 15)]
plt.plot([0],[0],marker='+',color = 'k',markersize = 20)
ax.legend(handles = legend_el, bbox_to_anchor=(0.5, 0.05),loc = 'lower right',frameon = False,facecolor= 'w', columnspacing = 10)
