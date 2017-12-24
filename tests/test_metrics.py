import unittest
from metrics import metrics as met
import pandas as pd
import numpy as np


class MetricsTests(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(MetricsTests, self).__init__(*args, **kwargs)
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

    def test_metrics(self):
        pass

    def test_auc(self):
        auc = self.metrics.auc(self.rawscore, self.newProg)
        actual_auc = 0.771
        self.assertAlmostEqual(actual_auc, auc)

    def test_bac(self):
        bac = self.metrics.bac(self.highrisk, self.newProg)
        actual_bac = 0.6353276
        self.assertAlmostEqual(actual_bac, bac)

    def test_mcc(self):
        mcc = self.metrics.mcc(self.highrisk, self.newProg)
        actual_mcc = 0.3833108
        self.assertAlmostEqual(actual_mcc, mcc)

    def test_f1(self):
        f1 = self.metrics.f1(self.highrisk, self.newProg)
        actual_f1 = 0.4444444
        self.assertAlmostEqual(actual_f1, f1)

    def test_timeroc(self):
        timeROC = self.metrics.timeROC(self.rawscore, self.PFStime, self.progression)
        # pd.DataFrame(index=[t=427     t=488     t=549     t=610     t=671])
        actual_timeROC = [0.7290861, 0.8157194, 0.7675568, 0.7460048, 0.7118679]
        print("timeROC", timeROC)
        print("actualTimeROC", actual_timeROC)
        for i in range(0, len(actual_timeROC)):
            self.assertAlmostEqual(actual_timeROC[i], timeROC[i])

    def test_iauc(self):
        iauc = self.metrics.integratedAUC(self.rawscore, self.PFStime, self.progression)
        actual_iAUC = 0.7241948
        self.assertAlmostEqual(actual_iAUC, iauc)

    def test_prauc(self):
        prauc = self.metrics.prAUC(self.highrisk, self.newProg)
        actual_prAUC = 0.6835078
        self.assertAlmostEqual(actual_prAUC, prauc)

