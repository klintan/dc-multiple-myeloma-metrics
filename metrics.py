from sklearn.metrics import auc, recall_score, matthews_corrcoef, f1_score, average_precision_score
from lifelines.utils import concordance_index
import numpy as np
import timeROC as timeROCPackage

class Calculate():
    def __init__(self):
        pass

    ## Calculate the Area Under the ROC Curve
    def auc(self, predicted, actual):
        return auc(actual, predicted)

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
        tempAUC = timeROCPackage(T=D_PFS, delta=D_PFS_FLAG, marker=predicted, cause=1, times=times)
        pass

    def integratedAUC(self, predicted, D_PFS, D_PFS_FLAG, times=30.5 * np.asarray([14, 16, 18, 20, 22])):
        print(times)
        pass

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
        pass
        # HR = rep(NA, len(rawscore))
        # HR[pfs_flag == 1 & PFStime < cutoff * 30.5] < - 1
        # HR[PFStime >= cutoff * 30.5] < - 0
        # newProg = HR
        # progression = pfs_flag

        # calculate = Calculate()
        #
        # auc = calculate.auc(rawscore, newProg)
        # bac = calculate.bac(highrisk, newProg)
        # mcc = calculate.mcc(highrisk, newProg)
        # f1 = calculate.f1(highrisk, newProg)
        # timeROC = calculate.timeROC(rawscore, PFStime, progression)
        # iAUC = calculate.integratedAUC(rawscore, PFStime, progression, times=30.5 * seq(12, 24, by=.25))
        #
        # # Remove null
        # rawscore_nona = rawscore[! is.na(HR)]
        # newProg_nona = newProg[! is.na(HR)]
        # prAUC = calculate.prAUC(rawscore_nona, newProg_nona)
        #
        # return (list(auc, bac, mcc, f1, timeROC, iAUC, prAUC))

    ##### simple wrapper to faciliate bootstrapping later
    def weightedMetrics(self, singleSubPredMat, PFStime, pfs_flag, study):
        pass
        # rawscore < - singleSubPredMat$predictionscore
        # highrisk < - as.numeric(as.logical(singleSubPredMat$highriskflag));
        # cutoff < - 18
        # HR < - rep(NA, length(rawscore))
        # HR[pfs_flag == 1 & PFStime < cutoff * 30.5] < - 1
        # HR[PFStime >= cutoff * 30.5] < - 0
        # newProg < - HR
        # progression < - pfs_flag
        # N < - table(study[progression == 1])
        #
        # auc < - c();
        # bac < - c();
        # mcc < - c();
        # f1 < - c();
        # timeROC < - c();
        # iAUC < - c();
        # prAUC < - c();
        #
        # for (s in names(N))
        # {
        #                     inds < - study == s
        # auc < - c(auc, calculate.auc(rawscore[inds], newProg[inds]))
        # bac < - c(bac, calculate.bac(highrisk[inds], newProg[inds]))
        # mcc < - c(mcc, calculate.mcc(highrisk[inds], newProg[inds]))
        # f1 < - c(f1, calculate.f1(highrisk[inds], newProg[inds]))
        # timeROC < - c(timeROC, calculate.timeROC(rawscore[inds], PFStime[inds], progression[inds]))
        # iAUC < - c(iAUC, calculate.integratedAUC(rawscore[inds], PFStime[inds], progression[inds], times =  30.5 * seq(12, 24, by = .25)))
        #
        # # Remove nulls for prAUC
        # rawscore_nona < - rawscore[inds][! is.na(HR[inds])]
        # newProg_nona < - newProg[inds][! is.na(HR[inds])]
        # prAUC < - c(prAUC, calculate.prAUC(rawscore_nona, newProg_nona))
        # }
        #
        # auc < - calculate.weightedAverage(auc, N = as.vector(N))
        # bac < - calculate.weightedAverage(bac, N = as.vector(N))
        # mcc < - calculate.weightedAverage(mcc, N = as.vector(N))
        # f1 < - calculate.weightedAverage(f1, N = as.vector(N))
        # timeROC < - calculate.weightedAverage(timeROC, N = as.vector(N))
        # iAUC < - calculate.weightedAverage(iAUC, N = as.vector(N))
        # prAUC < - calculate.weightedAverage(prAUC, N = as.vector(N))
        # return (list(auc, bac, mcc, f1, iAUC, prAUC))


