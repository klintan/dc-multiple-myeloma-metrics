import metrics
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from ipcw import IPCW
calculate = metrics.Calculate()
import numpy as np
# load ground truth test data

test = pd.read_csv('test_data.csv')
pred = test

# MinMaxScale the prediction scores, to get values similar to probability/confidence values (not necessary)
scaler = MinMaxScaler()
pred['prediction'] = pred['D_PFS'].values.argsort()
pred['prediction'] = scaler.fit_transform(pred['prediction'].values.reshape(-1, 1))

# IPCW examples
# T              : vector of observed failure times
# delta          : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
T = test['D_PFS'].values
delta = test['D_PFS_FLAG'].values

#IPCW_data = np.array((T, delta)).transpose()
df_ipcw_data = test[['D_PFS', 'D_PFS_FLAG']]
df_ipcw_data.columns = ['failure_time', 'status']
times = 30.5 * np.asarray([14, 16, 18, 20, 22])
# df_ipcw_data = pd.DataFrame(data=IPCW_data, index=IPCW_data[:, 0], columns=['failure_time', 'status'])

# remove duplicates using the general technique described by kaplan meier paper (perhaps we need to move this to
# somewhere else.
for i, row  in df_ipcw_data[df_ipcw_data['failure_time'].duplicated()].iterrows():
    if row['status'] == 0:
        df_ipcw_data.at[i, 'failure_time'] = row['failure_time']+1
    else:
        df_ipcw_data.at[i, 'failure_time'] = row['failure_time']-1

ipcw = IPCW(formula=None, data=df_ipcw_data, method="marginal", times=times, subjectTimes=T,
            what=["IPCW.times", "IPCW.subject.times"], subjectTimesLag=1)

weights = ipcw.marginal()
print(weights)
# predicted: continuous prediction score for high risk flag (patient progression within < 18 month.)
#calculate.timeROC(predicted=pred['D_PFS'], D_PFS=test['D_PFS'], D_PFS_FLAG=test['D_PFS_FLAG'])

# predict all metrics
#calculate.weightedMetrics(singleSubPredMat=pred, PFStime=test['D_PFS'], pfs_flag=test['D_PFS_FLAG'])
