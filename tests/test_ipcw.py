from metrics.ipcw.ipcw import IPCW
import numpy as np
import pandas as pd
import unittest

class IPCWTests(unittest.TestCase):
    def test_ipcw_marginal(self):
        # IPCW examples
        # T              : vector of observed failure times
        # delta          : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
        test = pd.read_csv('../sample-data/test_data.csv')

        T = test['D_PFS'].values
        delta = test['D_PFS_FLAG'].values

        IPCW_data = np.array((T, delta)).transpose()
        df_ipcw_data = test[['D_PFS', 'D_PFS_FLAG']]
        df_ipcw_data.columns = ['failure_time', 'status']
        times = 30.5 * np.asarray([14, 16, 18, 20, 22])
        ipcw = IPCW(formula=None, data=df_ipcw_data, method="marginal", times=times, subjectTimes=T,
                    what=["IPCW.times", "IPCW.subject.times"], subjectTimesLag=1)

        weights = ipcw.marginal()
        r_ipcw_weights =[0.839, 0.787, 0.760, 0.746, 0.746]
        print("weights", weights.values)
        print("actual weights", r_ipcw_weights)
        for i in range(0,len(r_ipcw_weights)):
            self.assertAlmostEqual(weights.values[i],r_ipcw_weights[i])

        #head() of the predicted IPCW for 5 subjects (rows), at the 1 requested times (columns):

        #[1] 0.839 0.787 0.760 0.746 0.746

        #head() of predicted IPCW at the individual subject times:

        #[1] 1.00 0.99 0.98 0.97 0.97 0.97

if __name__ == '__main__':
    unittest.main()