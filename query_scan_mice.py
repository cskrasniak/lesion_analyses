from oneibl.one import ONE
import numpy as np
import pandas as pd
from ibl_pipeline import subject, behavior
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
    data_DF['signedContrast'] = np.subtract(data_DF.contrastRight, data_DF.contrastLeft)
    data_DF['rt'] = np.subtract(data_DF.response_times, data_DF.goCueTimes)
    return data_DF


def DJ_fetch_DF(subjects, useDates):
    '''Fetches behavioral data (single trials) from dataJoint and returns it as a dataframe. You
    can restrict based on the subjects to use and the dates to retrieve data from. By default, it
    retrieves the subject uuid, session date, trial ID, time of response, choice, stim on time,
    signed contrast, and feedback type.
    subjects: list of strings, of subject names you wish to retrieve data for
    useDates: if a string of a single date, return all sessions including and after that date,
    otherwise a list of dates in format 'YYYY-MM-DD' that you wish to include'''
    DF_list = []

    for sub in subjects:
        print('grabbing data from ' + sub)
        subs = subject.Subject & 'subject_nickname = "{}"'.format(sub)

        if isinstance(useDates, str):
            trials = behavior.TrialSet.Trial & subs & ...
            'DATE(session_start_time) >= "{}"'.format(useDates)
            trials = trials * trials.proj(session_start_date='DATE(session_start_time)')
            trials = trials * trials.proj(
                signed_contrast='trial_stim_contrast_right - trial_stim_contrast_left')

            allSessions = pd.DataFrame(trials.fetch('subject_uuid', 'session_start_date',
                                                    'trial_id', 'trial_response_time',
                                                    'trial_response_choice',
                                                    'trial_stim_on_time', 'signed_contrast',
                                                    'trial_feedback_type', as_dict=True))
        elif isinstance(useDates, list):
            trials = behavior.TrialSet.Trial & subs & 'DATE(session_start_time) \
                >= "{}"'.format(useDates[0])
            trials = trials * trials.proj(
                signed_contrast='trial_stim_contrast_right - trial_stim_contrast_left')
            trials = trials * trials.proj(session_start_date='DATE(session_start_time)')
            allSessions = pd.DataFrame(trials.fetch('subject_uuid', 'session_start_date',
                                                    'trial_id', 'trial_response_time',
                                                    'trial_response_choice', 'trial_stim_on_time',
                                                    'signed_contrast', 'trial_feedback_type',
                                                    as_dict=True))

            useSessions = [ses.strftime('%Y-%m-%d') in useDates
                           for ses in allSessions.session_start_date]
            allSessions = allSessions[useSessions]

        allSessions['trial_response_choice'][allSessions['trial_response_choice'] == 'CW'] = 1
        allSessions['trial_response_choice'][allSessions['trial_response_choice'] == 'CCW'] = -1
        allSessions['trial_response_choice'][allSessions['trial_response_choice'] == 'No Go'] = 0

        allSessions['subject'] = sub
        DF_list.append(allSessions)
    return DF_list
