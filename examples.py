import metrics
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

calculate = metrics.Calculate()

# load ground truth test data

test = pd.read_csv('test_data.csv')
pred = test

# MinMaxScale the prediction scores, to get values similar to probability/confidence values (not necessary)
scaler = MinMaxScaler()
pred['prediction'] = pred['D_PFS'].values.argsort()
pred['prediction'] = scaler.fit_transform(pred['prediction'].values.reshape(-1, 1))

# predicted: continuous prediction score for high risk flag (patient progression within < 18 month.)
calculate.timeROC(predicted=pred['HR_FLAG'], D_PFS=test['D_PFS'], D_PFS_FLAG=test['D_PFS_FLAG'])
