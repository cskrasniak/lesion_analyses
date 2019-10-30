from oneibl.one import ONE
import numpy as np
import pandas as pd
one = ONE()


def ONE_query_scan_mice(subject, date_range):
    """ A function to get relevant ALYX data for sessions from a date range for a single subject
    subject = string of subject name, date_range = list of length 2 in format ['YYYY-MM-DD',
    'YYYY-MM-DD']. Returns: a pandas dataframe with all trials from sessions in that date range"""
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
