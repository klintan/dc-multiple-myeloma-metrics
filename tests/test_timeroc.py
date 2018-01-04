from metrics.timeroc.timeROC import timeROC
import pandas as pd
import numpy as np
from metrics import metrics as met
from metrics.timeroc.timeROC import timeROC, calculate_weights_cases_all, case_masking, matrix_transpose
import unittest
from rpy2 import robjects

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
        self.times = 30.5 * np.asarray([14, 16, 18, 20, 22])
        self.results = timeROC(T=self.PFStime, delta=self.progression, marker=self.rawscore, cause=1, times=self.times)

    def test_timeroc(self):
        # T              : vector of observed failure times
        # delta          : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
        # marker         : vector of marker values
        # cause          : the value that indicates the main event of interest
        # times          : vector of times you want to compute the time dependent AUC.
        timeROC = self.metrics.timeROC(self.rawscore, self.PFStime, self.progression)
        actual_timeROC = [0.7290861, 0.8157194, 0.7675568, 0.7460048, 0.7118679]
        print("timeROC", timeROC)
        print("actualTimeROC", actual_timeROC)
        for i in range(0, len(actual_timeROC)):
            self.assertAlmostEqual(actual_timeROC[i], timeROC[i])

    def test_tp(self):
        actual_tp = pd.read_csv("timeROC_tp.csv")
        tp = pd.DataFrame(self.results['TP'], columns=actual_tp.columns)
        print(actual_tp)
        print(tp)
        self.assertEqual(len(actual_tp.index), len(tp.index))

    def test_fp(self):
        actual_fp = pd.read_csv("timeROC_fp.csv")
        fp = pd.DataFrame(self.results['FP'], columns=actual_fp.columns)
        print(actual_fp)
        print(fp)
        self.assertEqual(len(actual_fp.index), len(fp.index))

    def test_survprob(self):
        actual_survprob = [0.7862705, 0.7365341, 0.7106746, 0.6838567, 0.6570388]
        survprob = self.results['survProb']
        self.assertTrue(all(np.isclose(survprob, actual_survprob)))

    def test_cumulative_incidence(self):
        actual_ci = [0.2019941, 0.2634659, 0.2893254, 0.3027343, 0.3429612]
        ci = self.results['CumulativeIncidence']
        self.assertTrue(all(np.isclose(ci, actual_ci)))

    def test_stats(self):
        actual_stats = pd.read_csv("timeROC_stats.csv")
        stats = self.results['Stats']
        print("actual", actual_stats)
        print(stats)

    def test_weights_cases_all(self):
        actual_weights_cases_all = [0.01041895, 0.01030928, 0.01400491, 0.01041895, 0.01156530, 0.01173538, 0.01191319,
                                    0.01248363, 0.01054299, 0.01228548, 0.07327658, 0.01030928, 0.01041895, 0.01054299,
                                    0.01248363, 0.01248363, 0.01400491, 0.01041895, 0.01609150, 0.01041895, 0.01041895,
                                    0.02630441, 0.01054299, 0.01094849, 0.01661058, 0.01269886, 0.01340895, 0.01041895,
                                    0.01124844, 0.01316064, 0.01517198, 0.02203807, 0.01400491, 0.01209647, 0.01400491,
                                    0.01437346, 0.01094849, 0.01173538, 0.01020408, 0.01156530, 0.01340895, 0.01400491,
                                    0.01779998, 0.01851198, 0.01292563, 0.01400491, 0.01716426, 0.03799526, 0.01476193,
                                    0.01228548, 0.01030928, 0.01156530, 0.01609150, 0.02442553, 0.02442553, 0.01041895,
                                    0.01109644, 0.09770210, 0.02012171, 0.01041895, 0.01010101, 0.01054299, 0.14655316,
                                    0.01340895, 0.29310631, 0.01248363, 0.01054299, 0.01561822, 0.01340895, 0.01080990,
                                    0.01370045, 0.01400491, 0.03419574, 0.03108703, 0.04885105, 0.01000000, 0.01094849,
                                    0.01269886, 0.01316064, 0.01340895, 0.01517198, 0.01779998, 0.01928331, 0.02313997,
                                    0.02849645, 0.01041895, 0.01140467, 0.01340895, 0.01716426, 0.02442553, 0.01094849,
                                    0.05862126, 0.04274467, 0.02442553, 0.02442553, 0.01716426, 0.02103634, 0.01340895,
                                    0.02313997, 0.04885105]
        T = self.PFStime.values
        n = len(T)
        delta = self.progression.values
        marker = self.rawscore.values
        order_marker = np.argsort(marker)

        Mat_data = pd.DataFrame(matrix_transpose(T, delta, marker), columns=["T", "delta", "marker"])
        Mat_data.sort_values("marker", ascending=False, inplace=True)

        weights_cases_all = calculate_weights_cases_all(Mat_data, order_marker, n)
        print("actual", actual_weights_cases_all)
        print(weights_cases_all)
        self.assertEqual(len(weights_cases_all), len(actual_weights_cases_all))
        self.assertTrue(all(np.isclose(weights_cases_all, actual_weights_cases_all)))

    def test_mat_data(self):
        T = self.PFStime.values
        delta = self.progression.values
        marker = self.rawscore.values
        Mat_data = pd.DataFrame(matrix_transpose(T, delta, marker), columns=["T", "delta", "marker"])
        Mat_data.sort_values("marker", ascending=False, inplace=True)

        actual_Mat_data = pd.read_csv("mat_data.csv")
        # print(actual_Mat_data.values)
        # print(Mat_data.values)

        # np.testing.assert_almost_equal(Mat_data['T'].values, actual_Mat_data['T'].values, decimal=0, err_msg='', verbose=True)
        # self.assertTrue(np.allclose(Mat_data['T'].values, actual_Mat_data['T'].values, atol=1e-02))
        # self.assertTrue(np.allclose(Mat_data['delta'].values, actual_Mat_data['delta'].values, atol=1e-05))
        # self.assertTrue(np.allclose(Mat_data['marker'].values, actual_Mat_data['marker'].values, atol=1e-05))
        for idx, val in enumerate(Mat_data['T'].values):
            diff = val - actual_Mat_data['T'].values[idx]
            print(diff)
            self.assertEqual(val, actual_Mat_data['T'].values[idx])

    def test_marker_sort(self):
        # I had a suspicion that the sorting in R and python differed and it makes a huge difference (if not only
        # because the other tests fail if its not correct)

        # timeROC.r line 68
        # original input marker values
        actual_marker = [28.360091, 3.505252, 3.319087, 28.360091, 5.156380, 1.000000, 51.563801,
                         9.471605, 35.625520, 8.414815, 4.421088, 15.978442, 77.102888, 85.909469,
                         10.000000, 50.683143, 9.647737, 64.773674, 15.282317, 51.563801, 48.921827,
                         24.000833, 6.287092, 7.005762, 25.938281, 28.844453, 6.653499, 47.160510,
                         78.864204, 40.115245, 6.125104, 26.907005, 11.407421, 36.594244, 61.251041,
                         29.813177, 51.563801, 64.773674, 6.477367, 50.683143, 82.386837, 28.360091,
                         5.068314, 2.250208, 5.156380, 3.130866, 29.813177, 64.773674, 5.156380,
                         63.012357, 28.360091, 5.156380, 6.477367, 15.282317, 27.875729, 28.360091,
                         60.554735, 47.250208, 87.670786, 51.563801, 51.563801, 40.469140, 29.813177,
                         5.068314, 2.069604, 13.344869, 51.563801, 49.187656, 5.596709, 3.054587,
                         28.360091, 6.125104, 12.693956, 24.969557, 5.156380, 33.688073, 3.659261,
                         26.907005, 5.948972, 5.068314, 5.420578, 5.260524, 16.251041, 43.637878,
                         1.017287, 17.219765, 5.156380, 9.382874, 36.594244, 5.156380, 13.344869,
                         2.074077, 5.156380, 5.156380, 48.921827, 5.156380, 4.005238, 28.844453,
                         36.592612, 35.625520]

        # timeROC.r line 124
        # order_T < - order(T)
        actual_order_T = [x - 1 for x in [5, 17, 42, 14, 20, 86, 29, 21, 68, 15, 59, 66, 43, 61, 34, 16, 58, 50,
                                          60, 8, 39, 31, 26, 77, 90, 28, 91, 89, 80, 51, 13, 83, 48, 98, 18, 47,
                                          35, 22, 37, 38, 7, 27, 45, 62, 55, 49, 99, 94, 10, 65, 30, 71, 23, 64,
                                          72, 36, 67, 1, 79, 41, 32, 76, 63, 74, 96, 100, 53, 12, 40, 84, 11, 44,
                                          78, 52, 56, 4, 93, 33, 92, 9, 75, 85, 70, 54, 19, 97, 46, 95, 87, 81,
                                          69, 25, 3, 6, 82, 2, 57, 73, 88, 24]]

        # ordered on T
        # timeROC.r line 127
        # marker <- marker[order_T]
        actual_marker_sorted_by_T = [5.156380, 9.647737, 28.360091, 85.909469, 51.563801, 17.219765, 78.864204,
                                     48.921827, 49.187656, 10.000000, 87.670786, 13.344869, 5.068314, 51.563801,
                                     36.594244, 50.683143, 47.250208, 63.012357, 51.563801, 9.471605, 6.477367,
                                     6.125104, 28.844453, 3.659261, 5.156380, 47.160510, 13.344869, 36.594244,
                                     5.068314, 28.360091, 77.102888, 16.251041, 64.773674, 28.844453, 64.773674,
                                     29.813177, 61.251041, 24.000833, 51.563801, 64.773674, 51.563801, 6.653499,
                                     5.156380, 40.469140, 27.875729, 5.156380, 36.592612, 5.156380, 8.414815,
                                     2.069604, 40.115245, 28.360091, 6.287092, 5.068314, 6.125104, 29.813177,
                                     51.563801, 28.360091, 5.948972, 82.386837, 26.907005, 33.688073, 29.813177,
                                     24.969557, 5.156380, 35.625520, 6.477367, 15.978442, 50.683143, 43.637878,
                                     4.421088, 2.250208, 26.907005, 5.156380, 28.360091, 28.360091, 5.156380,
                                     11.407421, 2.074077, 35.625520, 5.156380, 1.017287, 3.054587, 15.282317,
                                     15.282317, 4.005238, 3.130866, 48.921827, 5.156380, 5.420578, 5.596709,
                                     25.938281, 3.319087, 1.000000, 5.260524, 3.505252, 60.554735, 12.693956,
                                     9.382874, 7.005762]

        # timeROC.r line 151
        # order_marker <- order(- marker)
        # R has 1-indexing for arrays
        actual_order_marker_reversed = [x - 1 for x in
                                                   [11, 4, 60, 7, 31, 33, 35, 40, 18, 37, 97, 5, 14, 19, 39, 41, 57, 16,
                                                    69, 9, 8, 88, 17, 26, 70, 44, 51, 15, 28, 47, 66, 80, 62, 36, 56,
                                                    63, 23, 34, 3, 30, 52, 58, 75, 76, 45, 61, 73, 92, 64, 38, 6, 32,
                                                    68, 84, 85, 12, 27, 98, 78, 10, 2, 20, 99, 49, 100, 42, 21, 67, 53,
                                                    22, 55, 59, 91, 90, 95, 1, 25, 43, 46, 48, 65, 74, 77, 81, 89, 13,
                                                    29, 54, 71, 86, 24, 96, 93, 87, 83, 72, 79, 50, 82, 94]]

        # timeROC.r line 157
        # Mat_data <- cbind(T, delta, marker)[order_marker,]
        actual_mat_data_marker_sorted_by_order_marker = [87.670786, 85.909469, 82.386837, 78.864204, 77.102888,
                                                         64.773674, 64.773674, 64.773674, 63.012357, 61.251041,
                                                         60.554735, 51.563801, 51.563801, 51.563801, 51.563801,
                                                         51.563801, 51.563801, 50.683143, 50.683143, 49.187656,
                                                         48.921827, 48.921827, 47.250208, 47.160510, 43.637878,
                                                         40.469140, 40.115245, 36.594244, 36.594244, 36.592612,
                                                         35.625520, 35.625520, 33.688073, 29.813177, 29.813177,
                                                         29.813177, 28.844453, 28.844453, 28.360091, 28.360091,
                                                         28.360091, 28.360091, 28.360091, 28.360091, 27.875729,
                                                         26.907005, 26.907005, 25.938281, 24.969557, 24.000833,
                                                         17.219765, 16.251041, 15.978442, 15.282317, 15.282317,
                                                         13.344869, 13.344869, 12.693956, 11.407421, 10.000000,
                                                         9.647737, 9.471605, 9.382874, 8.414815, 7.005762, 6.653499,
                                                         6.477367, 6.477367, 6.287092, 6.125104, 6.125104, 5.948972,
                                                         5.596709, 5.420578, 5.260524, 5.156380, 5.156380, 5.156380,
                                                         5.156380, 5.156380, 5.156380, 5.156380, 5.156380, 5.156380,
                                                         5.156380, 5.068314, 5.068314, 5.068314, 4.421088, 4.005238,
                                                         3.659261, 3.505252, 3.319087, 3.130866, 3.054587, 2.250208,
                                                         2.074077, 2.069604, 1.017287, 1.000000]

        T = self.PFStime.values
        marker = self.rawscore.values

        order_T = np.argsort(T)
        marker_sorted_by_T = list(marker[order_T])

        marker_order = np.argsort(marker)
        marker_order_reversed = list(np.argsort(-marker))
        marker_sorted_by_marker = list(marker[marker_order])

        # check input to sort (redundant test, see input test)
        self.assertTrue(all(np.isclose(marker, actual_marker)))
        self.assertTrue(len(actual_order_T) - sum(np.isclose(order_T, actual_order_T)) < 5)

        # test the marker values sorted by T
        self.assertTrue(
            len(actual_marker_sorted_by_T) - sum(np.isclose(marker_sorted_by_T, actual_marker_sorted_by_T)) < 5)
        # check the sorted T values
        self.assertTrue(len(actual_order_T) - sum(np.isclose(order_T, actual_order_T)) < 5)

        self.assertTrue(len(actual_order_marker_reversed) - sum(np.isclose(marker_order_reversed, actual_order_marker_reversed)) < 5)

        #self.assertSequenceEqual(marker_order, actual_marker_order)

    def test_compare_rorder_pyargsort(self):
        rorder =  robjects.r['order']
        res = rorder(robjects.FloatVector(-np.array([28.360091, 3.505252, 3.319087, 28.360091, 5.156380, 1.000000, 51.563801,
                         9.471605, 35.625520, 8.414815, 4.421088, 15.978442, 77.102888, 85.909469,
                         10.000000, 50.683143, 9.647737, 64.773674, 15.282317, 51.563801, 48.921827,
                         24.000833, 6.287092, 7.005762, 25.938281, 28.844453, 6.653499, 47.160510,
                         78.864204, 40.115245, 6.125104, 26.907005, 11.407421, 36.594244, 61.251041,
                         29.813177, 51.563801, 64.773674, 6.477367, 50.683143, 82.386837, 28.360091,
                         5.068314, 2.250208, 5.156380, 3.130866, 29.813177, 64.773674, 5.156380,
                         63.012357, 28.360091, 5.156380, 6.477367, 15.282317, 27.875729, 28.360091,
                         60.554735, 47.250208, 87.670786, 51.563801, 51.563801, 40.469140, 29.813177,
                         5.068314, 2.069604, 13.344869, 51.563801, 49.187656, 5.596709, 3.054587,
                         28.360091, 6.125104, 12.693956, 24.969557, 5.156380, 33.688073, 3.659261,
                         26.907005, 5.948972, 5.068314, 5.420578, 5.260524, 16.251041, 43.637878,
                         1.017287, 17.219765, 5.156380, 9.382874, 36.594244, 5.156380, 13.344869,
                         2.074077, 5.156380, 5.156380, 48.921827, 5.156380, 4.005238, 28.844453,
                         36.592612, 35.625520])))
        print(res)



    def test_sort_order_T(self):
        T = self.PFStime.values
        order_T = np.argsort(T)

        # R has 1-indexing for arrays
        actual_order_T = [x - 1 for x in [5, 17, 42, 14, 20, 86, 29, 21, 68, 15, 59, 66, 43, 61, 34, 16, 58, 50,
                                          60, 8, 39, 31, 26, 77, 90, 28, 91, 89, 80, 51, 13, 83, 48, 98, 18, 47,
                                          35, 22, 37, 38, 7, 27, 45, 62, 55, 49, 99, 94, 10, 65, 30, 71, 23, 64,
                                          72, 36, 67, 1, 79, 41, 32, 76, 63, 74, 96, 100, 53, 12, 40, 84, 11, 44,
                                          78, 52, 56, 4, 93, 33, 92, 9, 75, 85, 70, 54, 19, 97, 46, 95, 87, 81,
                                          69, 25, 3, 6, 82, 2, 57, 73, 88, 24]]

        # the np.argsort does not sort exactly the same as R "order"-function, meaning we get some index swapping
        self.assertTrue(len(actual_order_T) - sum(np.isclose(order_T, actual_order_T)) < 5)

    def test_unique_marker_length(self):
        actual_marker = [28.360091, 3.505252, 3.319087, 28.360091, 5.156380, 1.000000, 51.563801,
                         9.471605, 35.625520, 8.414815, 4.421088, 15.978442, 77.102888, 85.909469,
                         10.000000, 50.683143, 9.647737, 64.773674, 15.282317, 51.563801, 48.921827,
                         24.000833, 6.287092, 7.005762, 25.938281, 28.844453, 6.653499, 47.160510,
                         78.864204, 40.115245, 6.125104, 26.907005, 11.407421, 36.594244, 61.251041,
                         29.813177, 51.563801, 64.773674, 6.477367, 50.683143, 82.386837, 28.360091,
                         5.068314, 2.250208, 5.156380, 3.130866, 29.813177, 64.773674, 5.156380,
                         63.012357, 28.360091, 5.156380, 6.477367, 15.282317, 27.875729, 28.360091,
                         60.554735, 47.250208, 87.670786, 51.563801, 51.563801, 40.469140, 29.813177,
                         5.068314, 2.069604, 13.344869, 51.563801, 49.187656, 5.596709, 3.054587,
                         28.360091, 6.125104, 12.693956, 24.969557, 5.156380, 33.688073, 3.659261,
                         26.907005, 5.948972, 5.068314, 5.420578, 5.260524, 16.251041, 43.637878,
                         1.017287, 17.219765, 5.156380, 9.382874, 36.594244, 5.156380, 13.344869,
                         2.074077, 5.156380, 5.156380, 48.921827, 5.156380, 4.005238, 28.844453,
                         36.592612, 35.625520]

        n_marker = len(set(actual_marker))
        actual_n_marker = 65

        self.assertEqual(n_marker, actual_n_marker)

    def test_masking(self):
        pass

    def test_cases(self):
        pass

    def test_control_1(self):
        pass

    def test_control_2(self):
        pass

    def test_timeroc_inputs(self):
        T = self.PFStime.values
        delta = self.progression.values
        marker = self.rawscore.values
        cause = 1
        times = self.times

        actual_T = [771.6500, 2407.4667, 1630.7333, 1128.5000, 10.0000, 1637.8500, 454.4500,
                    266.0000, 1165.1000, 610.0000, 1005.4833, 866.2000, 407.0000, 61.0000,
                    130.0000, 200.0000, 36.0000, 428.0000, 1314.5500, 61.0000, 100.6500,
                    441.0000, 701.5000, 2885.3000, 1602.2667, 292.0000, 459.0000, 348.0000,
                    93.0000, 669.9833, 274.0000, 782.8333, 1152.9000, 196.0000, 435.0000,
                    719.0000, 442.2500, 445.0000, 266.0000, 875.0000, 778.0000, 54.9000,
                    155.0000, 1013.6167, 521.5500, 1393.8500, 433.0000, 427.0000, 542.9000,
                    220.6167, 390.4000, 1110.2000, 854.0000, 1301.3333, 537.0000, 1125.4500,
                    2569.1167, 212.0000, 131.0000, 244.0000, 167.7500, 527.0000, 811.0000,
                    717.0000, 633.3833, 149.4500, 728.9500, 109.8000, 1572.0000, 1288.1167,
                    683.2000, 718.0000, 2626.0500, 812.0000, 1198.6500, 792.0000, 295.8500,
                    1101.0000, 774.0000, 386.0000, 1566.6833, 2006.9000, 421.0000, 982.0000,
                    1241.3500, 84.0000, 1512.8000, 2773.4667, 385.0000, 314.1500, 350.0000,
                    1163.0000, 1149.8500, 579.5000, 1469.0833, 820.4500, 1385.7167, 427.0000,
                    567.0000, 841.0000]
        actual_delta = [1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0,
                        1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0,
                        0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                        1, 0, 1, 1, 0, 0, 0]
        actual_marker = [28.360091, 3.505252, 3.319087, 28.360091, 5.156380, 1.000000, 51.563801,
                         9.471605, 35.625520, 8.414815, 4.421088, 15.978442, 77.102888, 85.909469,
                         10.000000, 50.683143, 9.647737, 64.773674, 15.282317, 51.563801, 48.921827,
                         24.000833, 6.287092, 7.005762, 25.938281, 28.844453, 6.653499, 47.160510,
                         78.864204, 40.115245, 6.125104, 26.907005, 11.407421, 36.594244, 61.251041,
                         29.813177, 51.563801, 64.773674, 6.477367, 50.683143, 82.386837, 28.360091,
                         5.068314, 2.250208, 5.156380, 3.130866, 29.813177, 64.773674, 5.156380,
                         63.012357, 28.360091, 5.156380, 6.477367, 15.282317, 27.875729, 28.360091,
                         60.554735, 47.250208, 87.670786, 51.563801, 51.563801, 40.469140, 29.813177,
                         5.068314, 2.069604, 13.344869, 51.563801, 49.187656, 5.596709, 3.054587,
                         28.360091, 6.125104, 12.693956, 24.969557, 5.156380, 33.688073, 3.659261,
                         26.907005, 5.948972, 5.068314, 5.420578, 5.260524, 16.251041, 43.637878,
                         1.017287, 17.219765, 5.156380, 9.382874, 36.594244, 5.156380, 13.344869,
                         2.074077, 5.156380, 5.156380, 48.921827, 5.156380, 4.005238, 28.844453,
                         36.592612, 35.625520]
        actual_cause = 1
        actual_times = [427, 488, 549, 610, 671]
        self.assertTrue(all(np.isclose(T, actual_T)))
        self.assertTrue(all(np.isclose(delta, actual_delta)))
        self.assertTrue(all(np.isclose(marker, actual_marker)))
        self.assertEqual(cause, actual_cause)
        self.assertTrue(all(np.isclose(times, actual_times)))
