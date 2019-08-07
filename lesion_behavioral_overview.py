import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from os.path import join
import seaborn as sns
import datajoint as dj
from ibl_pipeline import subject, acquisition, action, behavior, reference
from ibl_pipeline.analyses import behavior as behavior_analysis

# Settings
path = '/home/chris/DJ_figs/Behavior/'

# Query list of subjects
all_sub = subject.Subject * subject.SubjectLab & 'subject_birth_date > "2018-09-01"' & 'subject_line IS NULL OR subject_line="C57BL/6J"'
#all_sub = subject.Subject * subject.SubjectLab & 'subject_nickname = "ZM_1742"'
zadorSubs = pd.DataFrame(subject.Subject * subject.SubjectLab & 'lab_name = "zadorlab"')
subjects = pd.DataFrame(all_sub)

lesion_train_sub = pd.DataFrame(subject.Subject * subject.SubjectLab & 'subject_nickname = "CSK-les-006"' & 'subject_nickname = "CSK-les-007"' & 'subject_nickname = "CSK-les-008"')
lesion_then_train = ["CSK-les-006","CSK-les-007","CSK-les-008"]
train_then_lesion = ["CSK-IBL-003","CSK-IBL-005"]  
lesion_date = ["2019-04-07","2019-05-07"]      
learning = pd.DataFrame(columns=['mouse','lab','learned','date_learned','training_time','perf_easy','n_trials','training_trials','threshold','bias','reaction_time','lapse_low','lapse_high','lesion'])
for i, nickname in enumerate(subjects['subject_nickname']):
    print('Processing subject %s'%nickname)
    
    # Gather behavioral data for subject
    subj = subject.Subject & 'subject_nickname="%s"'%nickname
    behav = pd.DataFrame((behavior_analysis.BehavioralSummaryByDate * subject.Subject * subject.SubjectLab &
       'subject_nickname="%s"'%nickname).proj('session_date', 'performance_easy').fetch(as_dict=True, order_by='session_date'))
    behav_sesh = pd.DataFrame((behavior.TrialSet * subject.Subject * subject.SubjectLab &
       'subject_nickname="%s"'%nickname).proj('n_trials', 'n_correct_trials').fetch(as_dict=True, order_by='session_start_time'))
    rt = pd.DataFrame(((behavior_analysis.BehavioralSummaryByDate.ReactionTimeByDate * subject.Subject * subject.SubjectLab &
       'subject_nickname="%s"'%nickname)).proj('session_date', 'median_reaction_time').fetch(as_dict=True, order_by='session_date'))
    psych = pd.DataFrame(((behavior_analysis.BehavioralSummaryByDate.PsychResults * subject.Subject * subject.SubjectLab &
       'subject_nickname="%s"'%nickname)).proj('session_date', 'n_trials_stim','threshold','bias','lapse_low','lapse_high').fetch(as_dict=True, order_by='session_date'))
    
    # Find sessions to use
    
    first_trained_session = subj.aggr(behavior_analysis.SessionTrainingStatus &	'training_status="trained"', first_trained='DATE(min(session_start_time))')
    trainedSessions = pd.DataFrame((subj * behavior_analysis.SessionTrainingStatus * acquisition.Session & 'training_status = "trained"').proj('task_protocol', 'training_status').fetch(as_dict = True, order_by='session_start_time'))
    
        

    untrainable_session = subj.aggr(behavior_analysis.SessionTrainingStatus & 'training_status="untrainable"', first_trained='DATE(min(session_start_time))')
    ephys_session = subj.aggr(behavior_analysis.SessionTrainingStatus & 'training_status="ready for ephys"', first_trained='DATE(min(session_start_time))')
    if len(first_trained_session) == 0 & len(untrainable_session) == 0:
        learning.loc[i,'learned'] = 'in training'
        learning.loc[i,'training_time'] = len(behav)
    elif len(first_trained_session) == 0 & len(untrainable_session) == 1:
        learning.loc[i,'learned'] = 'untrainable'
        learning.loc[i,'training_time'] = len(behav)
    else:
        first_trained_session_dt = subj.aggr(behavior_analysis.SessionTrainingStatus &	'training_status="trained"', first_trained='min(session_start_time)').fetch1('first_trained')
        pretrainedSessions = pd.DataFrame((subj * behavior_analysis.SessionTrainingStatus * acquisition.Session & 'training_status = "training in progress"').proj('task_protocol', 'training_status').fetch(as_dict = True, order_by='session_start_time'))
        useSessions1 = pretrainedSessions.iloc[-3:None]
        trainedSessions['use']=['training' in task for task in trainedSessions['task_protocol']]
        useSessions2 = trainedSessions.loc[trainedSessions['use'] == True]
        del useSessions2['use']
        useSessions = useSessions2#useSessions1.append(useSessions2,ignore_index=True, sort =True)
        useSessions['session_date']= [i.date() for i in useSessions['session_start_time']]#takes all sessions where mouse is trained and performing training_CW and the three sessions that were used to define the mosue as trained
        first_trained_session_time = first_trained_session.fetch1('first_trained') 
        
        use_behav = behav[pd.Series([sesh in list(useSessions.session_date) for sesh in list(behav.session_date)])]#selecting only the days stipulated in useSessions
        use_psych = psych[pd.Series([sesh in list(useSessions.session_date) for sesh in list(psych.session_date)])]
        #use_behav_sesh = behav_sesh[pd.Series([sesh in list(useSessions.session_date) for sesh in list(psych.session_date)])]
        
        #trained_sessions = firstTrained_session
        learning.loc[i,'learned'] = 'trained'
        learning.loc[i,'date_learned'] = first_trained_session_time
        learning.loc[i,'training_time'] = sum(behav.session_date < first_trained_session_time) #get num sessions to trained
        learning.loc[i,'training_trials'] = sum(np.array(behav_sesh.loc[behav_sesh['session_start_time'] < first_trained_session_dt,['n_trials']]))#get num trials to trained

        learning.loc[i,'perf_easy'] = [[np.array(use_behav.performance_easy)]]
        learning.loc[i,'mean_perf_easy'] = np.mean(np.array(use_behav.performance_easy)) #get performance easy for all useSessions
        learning.loc[i,'n_useSessions'] = len(np.array(use_behav.performance_easy))
        learning.loc[i,'sd_perf_easy'] = np.std(np.array(use_behav.performance_easy))
        
        use_psych['n_trials'] = n_trials = [sum(s) for s in use_psych.n_trials_stim]
        learning.loc[i,'n_trials'] = [[np.array(use_psych.n_trials)]]
        learning.loc[i,'mean_n_trials'] = np.mean(np.array(use_psych.n_trials))
        learning.loc[i,'sd_n_trials'] = np.std(np.array(use_psych.n_trials))
        
        learning.loc[i,'threshold'] = [[np.array(use_psych.threshold)]]
        learning.loc[i,'mean_threshold'] = np.mean(np.array(use_psych.threshold))
        learning.loc[i,'sd_threshold'] = np.std(np.array(use_psych.threshold))
        
        learning.loc[i,'bias'] = [[np.array(use_psych.bias)]]
        learning.loc[i,'mean_bias'] = np.mean(np.array(use_psych.bias))
        learning.loc[i,'sd_bias'] = np.std(np.array(use_psych.bias))
        
        learning.loc[i,'lapse_low'] = [[np.array(use_psych.lapse_low)]]
        learning.loc[i,'mean_lapse_low'] = np.mean(np.array(use_psych.lapse_low))
        learning.loc[i,'sd_lapse_low'] = np.std(np.array(use_psych.lapse_low))
        
        learning.loc[i,'lapse_high'] = [[np.array(use_psych.lapse_high)]]
        learning.loc[i,'mean_lapse_high'] = np.mean(np.array(use_psych.lapse_high))
        learning.loc[i,'sd_lapse_high'] = np.std(np.array(use_psych.lapse_high))
        
        if nickname in lesion_then_train:
            learning.loc[i,'lesion'] = 1
        elif nickname in train_then_lesion:
            learning.loc[i,'lesion'] = 2
        else:
            learning.loc[i,'lesion'] = 0
        if sum(rt.session_date == first_trained_session_time) == 0:
            learning.loc[i,'reaction_time'] = float(rt.median_reaction_time[np.argmin(np.array(abs(rt.session_date - first_trained_session_time)))])*1000
        else:
            learning.loc[i,'reaction_time'] = float(rt.median_reaction_time[rt.session_date == first_trained_session_time])*1000
    if len(ephys_session) > 0:
        first_ephys_session_time = ephys_session.fetch1('first_trained')  
        learning.loc[i,'learned'] = 'ephys'
        learning.loc[i,'date_ephys'] = first_ephys_session_time
        learning.loc[i,'days_trained_ephys'] = sum((behav.session_date > first_trained_session_time) & (behav.session_date < first_ephys_session_time))
        
    # Add mouse info to dataframe
    learning.loc[i,'mouse'] = nickname
    learning.iloc[i]['lab'] = subjects.iloc[i]['lab_name']
    
# a weird bug makes a few of training_trials arrays, this is to fix that
for i in range(len(learning["training_trials"])):
    if type(learning.loc[i,"training_trials"]) == np.ndarray:
        learning.loc[i,'training_trials'] = learning.loc[i,'training_trials'][0]
# Select mice that learned
learned = learning[learning['learned'] == 'trained']
learned = learned.append(learning[learning['learned'] == 'ephys'])


# Convert to float
learned['training_time'] = learned['training_time'].astype(float)
learned['mean_perf_easy'] = learned['mean_perf_easy'].astype(float)
learned['mean_n_trials'] = learned['mean_n_trials'].astype(float)
learned['mean_threshold'] = learned['mean_threshold'].astype(float)
learned['mean_bias'] = learned['mean_bias'].astype(float)
learned['mean_lapse_low'] = learned['mean_lapse_low'].astype(float)
learned['mean_lapse_high'] = learned['mean_lapse_high'].astype(float)
learned['reaction_time'] = learned['reaction_time'].astype(float)

psychPlot = learned.loc[:,["mean_bias","mean_threshold",'mean_lapse_low','mean_lapse_high','lesion']]
plt.figure(1)
ax1 = plt.subplot(1,8,1)
#ax1 = sns.boxplot(data = psychPlot,y = "mean_bias", x='lesion')
ax1 = sns.stripplot(data = learned,y = "mean_bias", x='lesion', hue = 'mouse')
ax1.legend('')
sns.despine(right = True)

#plt.figure(2)
ax2 = plt.subplot(1,8,3)
#ax2 = sns.boxplot(data = psychPlot,y = "mean_threshold", x='lesion')
ax2 = sns.stripplot(data = learned,y = "mean_threshold", x='lesion', hue = 'mouse')
ax2.legend('')
sns.despine(right = True)

#plt.figure(3)
ax3 = plt.subplot(1,8,5)
#ax3 = sns.boxplot(data = psychPlot,y = "mean_lapse_low", x='lesion')
ax3 = sns.stripplot(data = learned,y = "mean_lapse_low", x='lesion', hue = 'mouse')
ax3.legend('')
sns.despine(right = True)

#plt.figure(4)
ax4 = plt.subplot(1,8,7)
#ax4 = sns.boxplot(data = psychPlot,y = "mean_lapse_high", x='lesion')
ax4 = sns.stripplot(data = learned,y = "mean_lapse_high", x='lesion', hue = 'mouse')
ax4.legend('')
sns.despine(right = True)

plt.figure(5)
ax5 = plt.subplot(1,8,1)
#ax5 = sns.boxplot(data = learned, y = 'training_time', x = 'lesion')
ax5 = sns.stripplot(data = learned, y = 'training_time', x = 'lesion', hue = 'mouse')
ax5.legend('')
sns.despine(right = True)

#plt.figure(6)
ax6 = plt.subplot(1,8,3)
#ax6 = sns.boxplot(data = learned, y = 'training_trials', x = 'lesion')
ax6 = sns.stripplot(data = learned, y = 'training_trials', x = 'lesion', hue = 'mouse')
ax6.legend('')
sns.despine(right = True)