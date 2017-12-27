from metrics.timeroc.timeROC import timeROC
import pandas as pd
import numpy as np
from metrics import metrics as met
from metrics.timeroc.timeROC import timeROC
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

    def test_tp(self):
        times = 30.5 * np.asarray([14, 16, 18, 20, 22])
        actual_tp = pd.read_csv("timeROC_tp.csv")
        results = timeROC(T=self.PFStime, delta=self.progression, marker=self.rawscore, cause=1, times=times)
        tp = pd.DataFrame(results['TP'], columns=actual_tp.columns)
        print(actual_tp)
        print(tp)
        self.assertEqual(len(actual_tp.index), len(tp.index))

    def test_fp(self):
        times = 30.5 * np.asarray([14, 16, 18, 20, 22])
        actual_fp = pd.read_csv("timeROC_fp.csv")
        results = timeROC(T=self.PFStime, delta=self.progression, marker=self.rawscore, cause=1, times=times)
        fp = pd.DataFrame(results['FP'], columns=actual_fp.columns)
        print(actual_fp)
        print(fp)
        self.assertEqual(len(actual_fp.index), len(fp.index))

    def test_survprob(self):
        times = 30.5 * np.asarray([14, 16, 18, 20, 22])
        actual_survprob = [0.7862705, 0.7365341, 0.7106746, 0.6838567, 0.6570388]
        results = timeROC(T=self.PFStime, delta=self.progression, marker=self.rawscore, cause=1, times=times)
        print("actual", actual_survprob)
        print(results['survProb'])

    def test_cumulative_incidence(self):
        times = 30.5 * np.asarray([14, 16, 18, 20, 22])
        actual_ci = [0.2019941, 0.2634659, 0.2893254, 0.3027343, 0.3429612]
        results = timeROC(T=self.PFStime, delta=self.progression, marker=self.rawscore, cause=1, times=times)
        ci = results['CumulativeIncidence']
        print("actual", actual_ci)
        print(ci)

    def test_stats(self):
        times = 30.5 * np.asarray([14, 16, 18, 20, 22])
        actual_stats = pd.read_csv("timeROC_stats.csv")
        results = timeROC(T=self.PFStime, delta=self.progression, marker=self.rawscore, cause=1, times=times)
        stats = results['Stats']
        print("actual", actual_stats)
        print(stats)
