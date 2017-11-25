## Calculate the Area Under the ROC Curve
calculate.auc <- function(predicted,actual) {
  suppressPackageStartupMessages(library("pROC"))
  pROC::auc(response=actual, predictor=predicted)
}

## Calculate the balanced accuracy [ i.e., ( sensitivity + specificity ) / 2 ]
calculate.bac <- function(predicted, actual) {
  suppressPackageStartupMessages(library("caret"))  
  cm <- confusionMatrix(data=predicted, reference=actual, positive = "1")
  as.numeric(cm$byClass["Balanced Accuracy"])
}

## Calculate Matthew's Correlation Coefficient
calculate.mcc <- function(predicted,actual) {
  suppressPackageStartupMessages(library("caret"))
  cm <- as.table(confusionMatrix(data= predicted,reference =actual), positive = "1")
  ((cm[1,1]*cm[2,2])-(cm[1,2]*cm[2,1]))/(sqrt(cm[1,1]+cm[1,2])*sqrt(cm[1,1]+cm[2,1])*sqrt(cm[2,2]+cm[1,2])*sqrt(cm[2,2]+cm[2,1]))
}

## Calculate F1 score
calculate.f1 <- function(predicted,actual) {
  suppressPackageStartupMessages(library("caret"))
  cm <- confusionMatrix(data=predicted,reference = actual, positive = "1")
  as.numeric(cm$byClass["F1"])
}

## Calculate TimeROC
calculate.timeROC <- function(predicted, D_PFS, D_PFS_FLAG, times = 30.5 * c(14, 16, 18, 20, 22)) {
  suppressPackageStartupMessages(library("timeROC"))
  suppressPackageStartupMessages(library("risksetROC"))
  tempAUC <- timeROC(T=D_PFS, delta=D_PFS_FLAG, marker=predicted,cause=1, times=times)
  iaucs <- IntegrateAUC(tempAUC$AUC, tempAUC$times, tempAUC$survProb,tmax = max(tempAUC$times))
  tAUCs <- tempAUC$AUC
}

## Calculate IntegratedAUC
calculate.integratedAUC <- function(predicted, D_PFS, D_PFS_FLAG, times = 30.5 * c(14, 16, 18, 20, 22)) {
  suppressPackageStartupMessages(library("timeROC"))
  suppressPackageStartupMessages(library("risksetROC"))
  tempAUC <- timeROC(T=D_PFS, delta=D_PFS_FLAG, marker=predicted,cause=1, times=times)
  iaucs <- IntegrateAUC(tempAUC$AUC, tempAUC$times, tempAUC$survProb,tmax = max(tempAUC$times))
  return(iaucs)
}

## Calculate precision recall based AUC
calculate.prAUC <- function(predicted, actual) {
  suppressPackageStartupMessages(library("PRROC"))
  prAUC <- pr.curve(scores.class0 = predicted, weights.class0 = actual)$auc.integral
}

## Calculate concordance Index ### package suport issues in ec2 instances. Should fix but low priority
calculate.concordanceIndex <- function(predicted, D_PFS, D_PFS_FLAG) {
  suppressPackageStartupMessages(library("survcomp"))
  cIndex <- survcomp::concordance.index(x = predicted, surv.time = D_PFS, surv.event  = D_PFS_FLAG)
}

## Wrapper for calculating all metrics: 
# rawscore    = continuous prediciton score from participant
# highrisk    = binary predicted score from participant
# newProg     = flag for true PFS_FLAG and PFS < 18mo
# PFStime     = true PFS time
# progression = true PFS_FLAGs

calculate.metrics <- function(rawscore, highrisk, PFStime, pfs_flag) { 

  cutoff                                    <- 18
  HR                                        <- rep(NA, length(rawscore))
  HR[pfs_flag == 1 & PFStime < cutoff*30.5] <- 1 
  HR[PFStime >= cutoff*30.5]                <- 0
  newProg                                   <- HR
  progression                               <- pfs_flag
  
  auc     <- calculate.auc(rawscore, newProg)
  bac     <- calculate.bac(highrisk, newProg)
  mcc     <- calculate.mcc(highrisk, newProg)
  f1      <- calculate.f1(highrisk, newProg)
  timeROC <- calculate.timeROC(rawscore, PFStime, progression) 
  iAUC    <- calculate.integratedAUC(rawscore, PFStime, progression, times =  30.5 * seq(12,24, by = .25))

  #Remove null
  rawscore_nona <- rawscore[!is.na(HR)]
  newProg_nona  <- newProg[!is.na(HR)]
  prAUC         <- calculate.prAUC(rawscore_nona, newProg_nona)
  
  return(list(auc,bac,mcc,f1,timeROC,iAUC,prAUC))
}

## function to calculate the weighted average, 
#  takes in array of one metric for each validation study (exampel auc) 
#  and teh N for each of the studies
#

calculate.weightedAverage <- function(metric, N)
{
  wAve   <- sum(metric*N)/sum(N)
  return(wAve)
}

##### simple wrapper to faciliate bootstrapping later
calculate.weightedMetrics <- function(singleSubPredMat, PFStime, pfs_flag, study) { 

  rawscore                                  <- singleSubPredMat$predictionscore 
  highrisk                                  <- as.numeric(as.logical(singleSubPredMat$highriskflag));
  cutoff                                    <- 18
  HR                                        <- rep(NA, length(rawscore))
  HR[pfs_flag == 1 & PFStime < cutoff*30.5] <- 1 
  HR[PFStime >= cutoff*30.5]                <- 0
  newProg                                   <- HR
  progression                               <- pfs_flag
  N                                         <- table(study[progression==1])
  
  auc <- c(); bac <- c(); mcc <- c(); f1 <- c(); timeROC <- c(); iAUC <- c(); prAUC <- c();
  
  for(s in names(N))
  {
    inds    <- study == s 
    auc     <- c(auc,calculate.auc(rawscore[inds], newProg[inds]))
    bac     <- c(bac,calculate.bac(highrisk[inds], newProg[inds]))
    mcc     <- c(mcc,calculate.mcc(highrisk[inds], newProg[inds]))
    f1      <- c(f1,calculate.f1(highrisk[inds], newProg[inds]))
    timeROC <- c(timeROC,calculate.timeROC(rawscore[inds], PFStime[inds], progression[inds]))
    iAUC    <- c(iAUC,calculate.integratedAUC(rawscore[inds], PFStime[inds], progression[inds], times =  30.5 * seq(12,24, by = .25)))

    #Remove nulls for prAUC
    rawscore_nona <- rawscore[inds][!is.na(HR[inds])]
    newProg_nona  <- newProg[inds][!is.na(HR[inds])]
    prAUC         <- c(prAUC,calculate.prAUC(rawscore_nona, newProg_nona))
  }

  auc     <- calculate.weightedAverage(auc,N = as.vector(N))
  bac     <- calculate.weightedAverage(bac,N = as.vector(N))
  mcc     <- calculate.weightedAverage(mcc,N = as.vector(N))
  f1      <- calculate.weightedAverage(f1,N = as.vector(N))
  timeROC <- calculate.weightedAverage(timeROC,N = as.vector(N))
  iAUC    <- calculate.weightedAverage(iAUC,N = as.vector(N))
  prAUC   <- calculate.weightedAverage(prAUC,N = as.vector(N))
  return(list(auc,bac,mcc,f1,iAUC,prAUC))
}







