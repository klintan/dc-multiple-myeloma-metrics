from metrics.timeroc.timeROC import timeROC
import pandas as pd
import numpy as np
from metrics import metrics as met
import unittest

class TimeROCTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(TimeROCTests, self).__init__(*args, **kwargs)
        test_data = pd.read_csv('../sample-data/test_data.csv')
        pred_data = pd.read_csv('../sample-data/predictions.csv')
        self.metrics = met.Calculate()

        self.PFStime = test_data['D_PFS']
        self.rawscore = pred_data['predictionscore']
        self.highrisk = pred_data['highriskflag']
        cutoff = 18
        self.HR = np.where((test_data['D_PFS'] < cutoff * 30.5) & (test_data['D_PFS_FLAG'] == 1), 1, 0)
        self.newProg = self.HR
        self.progression = test_data['D_PFS_FLAG']

    def test_timeroc(self):
        # T              : vector of observed failure times
        # delta          : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
        # marker         : vector of marker values
        # cause          : the value that indicates the main event of interest
        # times          : vector of times you want to compute the time dependent AUC.
        timeROC = self.metrics.timeROC(self.rawscore, self.PFStime, self.progression)
        # pd.DataFrame(index=[t=427     t=488     t=549     t=610     t=671])
        actual_timeROC = [0.7290861, 0.8157194, 0.7675568, 0.7460048, 0.7118679]
        print("timeROC", timeROC)
        print("actualTimeROC", actual_timeROC)
        for i in range(0, len(actual_timeROC)):
            self.assertAlmostEqual(actual_timeROC[i], timeROC[i])
