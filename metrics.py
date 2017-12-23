from sklearn.metrics import roc_auc_score, recall_score, matthews_corrcoef, f1_score, average_precision_score
from lifelines.utils import concordance_index
import numpy as np
from timeROC import timeROC
from integrateAUC import IntegrateAUC

class Calculate():
    def __init__(self):
        pass

    ## Calculate the Area Under the ROC Curve
    def auc(self, predicted, actual):
        return roc_auc_score(actual, predicted)

    ## Calculate the balanced accuracy [ i.e., ( sensitivity + specificity ) / 2 ]
    def bac(self, predicted, actual):
        # https://github.com/rhiever/tpot/issues/108
        # https://github.com/scikit-learn/scikit-learn/issues/6747
        return recall_score(actual, predicted, average='macro')

    ## Calculate Matthew's Correlation Coefficient
    def mcc(self, predicted, actual):
        return matthews_corrcoef(actual, predicted)

    ## Calculate F1 score
    def f1(self, predicted, actual):
        return f1_score(actual, predicted)

    ## Calculate precision recall based AUC
    def prAUC(self, predicted, actual):
        return average_precision_score(actual, predicted)

    ## Calculate TimeROC
    def timeROC(self, predicted, D_PFS, D_PFS_FLAG, times=30.5 * np.asarray([14, 16, 18, 20, 22])):
        tempAUC = timeROC(T=D_PFS, delta=D_PFS_FLAG, marker=predicted, cause=1, times=times)
        return tempAUC['AUC']

    def integratedAUC(self, predicted, D_PFS, D_PFS_FLAG, times=30.5 * np.asarray([14, 16, 18, 20, 22])):
        tempAUC = timeROC(T=D_PFS, delta=D_PFS_FLAG, marker=predicted, cause=1, times=times)
        iaucs = IntegrateAUC(tempAUC['AUC'], tempAUC['times'], tempAUC['survProb'], max(tempAUC['times']))
        return iaucs

    ## Calculate concordance Index
    def concordanceIndex(self, predicted, D_PFS, D_PFS_FLAG):
        return concordance_index(event_times=D_PFS, predicted_event_times=predicted, event_observed=D_PFS_FLAG)

    def weightedAverage(self, metric, N):
        ## function to calculate the weighted average,
        #  takes in array of one metric for each validation study (example auc)
        #  and the N for each of the studies
        return sum(metric * N) / sum(N)

    ## Wrapper for calculating all metrics:
    # rawscore    = continuous prediciton score from participant
    # highrisk    = binary predicted score from participant
    # newProg     = flag for true PFS_FLAG and PFS < 18mo
    # PFStime     = true PFS time
    # progression = true PFS_FLAGs
    def metrics(self, rawscore, highrisk, PFStime, pfs_flag):
        cutoff = 18
        HR1 = np.where((PFStime < cutoff * 30.5) & (pfs_flag == 1), 1, 0)
        HR2= np.where(PFStime >= cutoff * 30.5, 0, 1)
        assertCountEqual(HR1==HR2)
        newProg = HR
        progression = pfs_flag


        auc = self.auc(rawscore, newProg)
        bac = self.bac(highrisk, newProg)
        mcc = self.mcc(highrisk, newProg)
        f1 = self.f1(highrisk, newProg)
        timeROC = self.timeROC(rawscore, PFStime, progression)
        iAUC = self.integratedAUC(rawscore, PFStime, progression, times=30.5 * np.arange(12, 24, .25))

        # Remove null
        #rawscore_nona = rawscore[! is.na(HR)]
        #newProg_nona = newProg[! is.na(HR)]
        #prAUC = calculate.prAUC(rawscore_nona, newProg_nona)
        prAUC = 0

        return [auc, bac, mcc, f1, timeROC, iAUC, prAUC]

    ##### simple wrapper to faciliate bootstrapping later
    def weightedMetrics(self, singleSubPredMat, PFStime, pfs_flag, study=None):
        """
        :param singleSubPredMat: Full prediction matrix with columns: study,patient,predictionscore,highriskflag
        :param PFStime: actual time to failure
        :param pfs_flag: actual observed event flag
        :param study: the study str
        :return:
        """

        rawscore = singleSubPredMat['predictionscore']
        highrisk = singleSubPredMat['highriskflag']
        cutoff = 18
        HR = np.where((PFStime < cutoff * 30.5) & (pfs_flag == 1), 1, 0)
        # HR2 = np.where(PFStime >= cutoff * 30.5, 0, 1)
        newProg = HR
        progression = pfs_flag
        N = study[progression == 1]


        for s in singleSubPredMat['study'].unique():
            inds = singleSubPredMat['study'] == s
            auc = (self.auc(rawscore[inds], newProg[inds]))
            bac = self.bac(highrisk[inds], newProg[inds])
            mcc = self.mcc(highrisk[inds], newProg[inds])
            f1 = self.f1(highrisk[inds], newProg[inds])
            timeROC = self.timeROC(rawscore[inds], PFStime[inds], progression[inds])
            iAUC = self.integratedAUC(rawscore[inds], PFStime[inds], progression[inds], times =  30.5 * np.arange(12, 24, .25))

            # Remove nulls for prAUC
            #rawscore_nona = rawscore[inds][! is.na(HR[inds])]
            #newProg_nona =  newProg[inds][! is.na(HR[inds])]
            #prAUC = self.prAUC(rawscore_nona, newProg_nona)
            prAUC = []
            print("study: "+s, [auc, bac, mcc, f1, iAUC])

        #auc = self.weightedAverage(auc, N)
        #bac = self.weightedAverage(bac, N)
        #mcc = self.weightedAverage(mcc, N)
        #f1  = self.weightedAverage(f1, N)
        #timeROC = self.weightedAverage(timeROC, N)
        #iAUC = self.weightedAverage(iAUC, N)
        #prAUC = self.weightedAverage(prAUC, N)
        return [auc, bac, mcc, f1, iAUC, prAUC]


