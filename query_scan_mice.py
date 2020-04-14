from oneibl.one import ONE
import numpy as np
import pandas as pd
from ibl_pipeline import subject, behavior
import os
import platform
one = ONE()


def ONE_query_scan_mice(subject, date_range):
    """Get relevant ALYX data for sessions from a date range for a single subject
    subject = string of subject name
    date_range = list of length 2 in format ['YYYY-MM-DD','YYYY-MM-DD'].
    Returns: a pandas dataframe with all trials from sessions in that date range"""
    eid4 = one.search(subject=subject, date_range=date_range)
    dataset_types = ['trials.contrastRight', 'trials.feedbackType',
                     'trials.response_times', 'trials.contrastLeft',
                     'trials.choice', 'trials.goCueTrigger_times']
    data_DF = pd.DataFrame(columns=['contrastRight', 'feedback', 'response_times', 'contrastLeft',
                                    'choice', 'goCueTimes'])

    # loop over sessions
    for ses in range(len(eid4)):
        scan4 = one.load(eid4[ses], dataset_types=dataset_types, dclass_output=True)
        if scan4.dataset_id:
            DF = pd.DataFrame(np.transpose(np.array(scan4.data)),
                              columns=['contrastRight', 'feedback', 'response_times',
                                       'contrastLeft', 'choice', 'goCueTimes'])
            nans = np.isnan(DF)
            DF[nans] = 0  # make all nans in contrasts 0 so I can subtract them
            DF['EID'] = [eid4[ses]] * np.size(scan4.data, 1)
            data_DF = data_DF.append(DF)
    data_DF['signedContrast'] = np.subtract(data_DF.contrastLeft, data_DF.contrastRight)
    data_DF['rt'] = np.subtract(data_DF.response_times, data_DF.goCueTimes)
    return data_DF


def DJ_fetch_DF(subjects, useDates):
    '''Fetches behavioral data (single trials) from dataJoint and returns it as a dataframe. You
    can restrict based on the subjects to use and the dates to retrieve data from. By default, it
    retrieves the subject uuid, session date, trial ID, time of response, choice, stim on time,
    signed contrast, and feedback type.
    subjects: list of strings, of subject names you wish to retrieve data for
    useDates: if a string of a single date, return all sessions including and after that date,
    otherwise a list of dates in format 'YYYY-MM-DD' that you wish to include. So if you want
    data from a single day, useDates = ['YYYY-MM-DD']'''
    DF_list = []

    for sub in subjects:
       
        subs = subject.Subject & 'subject_nickname = "{}"'.format(sub)

        if isinstance(useDates, str):
            print('grabbing data from ' + sub + ' for ' + useDates)
            trials = behavior.TrialSet.Trial & subs & 'DATE(session_start_time) >= "{}"'.format(useDates)
            trials = trials * trials.proj(session_start_date='DATE(session_start_time)')
            trials = trials * trials.proj(signed_contrast='trial_stim_contrast_right - trial_stim_contrast_left')

            allSessions = pd.DataFrame(trials.fetch('subject_uuid', 'session_start_date',
                                                    'trial_id', 'trial_response_time',
                                                    'trial_response_choice',
                                                    'trial_go_cue_trigger_time', 'signed_contrast',
                                                    'trial_feedback_type', as_dict=True))
        elif isinstance(useDates, list):
            print('grabbing data from ' + sub + ' for ' + useDates[0])
            trials = behavior.TrialSet.Trial & subs & 'DATE(session_start_time) \
                >= "{}"'.format(useDates[0])
            trials = trials * trials.proj(
                signed_contrast='trial_stim_contrast_right - trial_stim_contrast_left')
            trials = trials * trials.proj(session_start_date='DATE(session_start_time)')
            allSessions = pd.DataFrame(trials.fetch('subject_uuid', 'session_start_date',
                                                    'trial_id', 'trial_response_time',
                                                    'trial_response_choice',
                                                    'trial_go_cue_trigger_time', 'signed_contrast',
                                                    'trial_feedback_type', as_dict=True))

            useSessions = [ses.strftime('%Y-%m-%d') in useDates
                           for ses in allSessions['session_start_date']]
            allSessions = allSessions[useSessions]

        #allSessions['trial_response_choice'][allSessions['trial_response_choice'] == 'CW'] = 1
        #allSessions['trial_response_choice'][allSessions['trial_response_choice'] == 'CCW'] = -1
        #allSessions['trial_response_choice'][allSessions['trial_response_choice'] == 'No Go'] = 0

        allSessions['subject'] = sub
        DF_list.append(allSessions)
    return DF_list


def align_laser2behavior(subjects):
    '''Takes the behavior from scanningBiasedChoiceWorld and aligns it to the laser position data
    that has been changed from a .m file to .npy and added to the subjects folder by using the
    laser2npy.m function. this takes in a list of subjects as an argument and outputs a list of
    dataframes that have simple behavior data as well as laser positions added to it'''
    if platform.system() == 'Darwin':
        dataPath = '/Users/ckrasnia/Desktop/Zador_Lab/scanData/Subjects'
    else:
        dataPath = "F:\Subjects"
    allData = []
   #  print('getting laser data from {}'.format(dataPath))
    for sub in subjects:  # loop over subjects
        os.chdir(dataPath)
        days = os.listdir(sub)
        training = []

        for day in days: # loop over days
            os.chdir(dataPath)
            os.chdir(sub)
            runs = os.listdir(day)
            behav = DJ_fetch_DF([sub], [day])
            behav = behav[0]

            for run in runs:
                print('Run #{}'.format(run))
                os.chdir(os.path.join(dataPath, sub, day, run))
                if len(os.listdir()) >= 1:  # if there is laser data in the folder, load it
                    ld = np.load("laserData")
# if the laser data is one longer than the behavior, remove the last laser and align to the start
                if len(ld) - 1 == np.size(behav, 0):
                    if np.shape(ld)[1] == 1:
                        behav['laserOn'] = ld[:-1, 0]
                        training.append(behav)
                    elif ld[0, 2].any():
                        behav['laserPosX'] = ld[:-1, 0]
                        behav['laserPosY'] = ld[:-1, 1]
                        training.append(behav)
                    else:
                        print('skipping session {} {}, marked as "laser off"'.format(day, run))

                elif len(ld) == np.size(behav, 0) - 1:
                    if np.shape(ld)[1] == 1:
                        behav.drop(behav.tail(1).index, inplace=True)
                        behav['laserOn'] = ld[:, 0]
                        training.append(behav)
                    elif ld[0, 2].any():             
                        behav.drop(behav.tail(1).index, inplace=True)
                        behav['laserPosX'] = ld[:, 0]
                        behav['laserPosY'] = ld[:, 1]
                        training.append(behav)
                    else:
                        print('skipping session {} {}, marked as "laser off"'.format(day, run))
                        
                elif len(ld) - 2 == np.size(behav, 0):
                    if np.shape(ld)[1] == 1:
                        behav['laserOn'] = ld[:-2, 0]
                        training.append(behav)
                    elif ld[0, 2].any():
                        behav['laserPosX'] = ld[:-1, 0]
                        behav['laserPosY'] = ld[:-1, 1]
                        training.append(behav)
                    else:
                        print('skipping session {} {}, marked as "laser off"'.format(day, run))

                elif len(ld) == np.size(behav, 0):
                    if np.shape(ld)[1] == 1:
                        behav['laserOn'] = ld[:, 0]
                        training.append(behav)
                    elif ld[0, 2].any():  # if marked 'laser on'
                        behav['laserPosX'] = ld[:, 0]
                        behav['laserPosY'] = ld[:, 1]
                        training.append(behav)
                    elif not ld[0, 2].all():
                        print('skipping session {} {}, marked as "laser off"'.format(day, run))
                    
                else:  # if they aren't matching or one off by the laser longer, I can't trust it
                    print('cannot align, for session {} laser #: {} trials #: {}'.format(day,
                          len(ld), np.size(behav, 0)))
              
        data = pd.concat(training, ignore_index=True)
        allData.append(data)
    return allData
