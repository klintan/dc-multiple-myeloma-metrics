import pandas as pd

#tp = pd.read_csv("timeROC_tp.csv", sep=" ")
#print(tp.head())
#fp = pd.read_csv("timeROC_fp.csv", sep=" ")
#fp = fp[['t=427', 't=488', 't=549', 't=610','t=671']].reset_index(drop=True)
#print(fp.head())
#fp.to_csv("timeROC_fp_new.csv", index=False)
#tp.to_csv("timeROC_tp_new.csv", index=False)

stats = pd.read_csv("timeROC_stats.csv", sep=" ")
stats.to_csv("timeROC_stats.csv", index=False)
print(stats)