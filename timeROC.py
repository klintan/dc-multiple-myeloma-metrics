import numpy as np
import functools

def reduce_concat(x, sep=""):
    return functools.reduce(lambda x, y: str(x) + sep + str(y), x)

def paste(*lists, sep=" ", collapse=None):
    result = map(lambda x: reduce_concat(x, sep=sep), zip(*lists))
    if collapse is not None:
        return reduce_concat(result, sep=collapse)
    return list(result)

def timeROC(T, delta, marker, cause, times, other_markers=None, weighting="marginal", ROC=True, iid=False):
    # https://github.com/cran/timeROC/blob/master/R/timeROC_3.R

    # T              : vector of observed failure times
    # delta          : vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
    # marker         : vector of marker values
    # other_markers  : (default is NULL, should be a matrix) other markers that can be associated with the censoring mechanism
    # cause          : the value that indicates the main event of interest
    # weighting      : (default is "marginal") weighting technique for IPCW : "marginal" for Kaplan-Meier, "cox" for proportional hazards cox model, "aalen" for additive aalen model
    # times          : vector of times you want to compute the time dependent AUC.
    # ROC            : if TRUE, then save True Positive fraction (Sensitivity) and False Positive fraction (1-Specificity)
    #                  for all time in vetor times
    # iid            : TRUE or FALSE, indicates if we want to compute the iid representation of the AUC estimator

    # check some inputs
    if len(delta) != len(T) or len(marker) != len(T) or len(marker) != len(delta):
        raise ("lengths of vector T, delta and marker have to be equal\n")

    if not times:
        raise ("Choose at least one time for computing the time-dependent AUC\n")

    if not weighting in ["marginal", "cox", "aalen"]:
        raise ("the weighting argument must be marginal (default), cox or aalen.\n")

    if weighting in ["cox", "aalen"] and isinstance(other_markers, (np.ndarray, np.generic)):
        raise ("argument other_markers must be numpy array/matrix\n")

    if weighting in ["cox", "aalen"] and other_markers:
        if np.shape(other_markers)[0] == len(marker):
            raise("lengths of vector T, delta, marker and number of rows of other_markers have to be equal\n")

    # check if there are missing values, and delete rows with missing values
    # not implemented yet

    # Initialize variables
    n = len(T)
    #n_marker = len(unique(marker)) # not sure what this is suppose to do
    n_times = len(times)
    if n_times == 1:
        times = [0] + times
        n_times = 2  # trick to use ipcw.cox() even if there is only one time

    times = sorted(times)
    times_names = paste("t=", times, sep="")

    # Outputs
    AUC_1 = np.zeros(n_times) # I don't think these are necessary but will stay for now
    AUC_2 = np.zeros(n_times) # I don't think these are necessary but will stay for now
    CumInci = np.zeros(n_times)
    surv = np.zeros(n_times)

    #names(AUC_1) < -times_names
    #names(AUC_2) < -times_names
    #names(CumInci) < -times_names
    #names(surv) < -times_names

    Stats < -matrix(NA, nrow=n_times, ncol=4)
    colnames(Stats) < -c("Cases", "survivor at t", "Other events at t", "Censored at t")
    rownames(Stats) < -times_names
    # }}}
    # {{{  computation of weights (1/2)
    # we need to order to use the pec::ipcw() fonction
    order_T < -order(T)
    T < - T[order_T]
    delta < - delta[order_T]
    marker < - marker[order_T]
    # use ipcw function from pec package
    if (weighting == "marginal"){
    weights < - pec::
        ipcw(Surv(failure_time, status)
    ~1, data = data.frame(failure_time=T, status=as.numeric(
        delta != 0)), method = "marginal", times = times, subjectTimes = T, subjectTimesLag = 1)
    }
    if (weighting == "cox"){
    if (missing(other_markers)){marker_censoring < -marker}
    other_markers < -other_markers[order_T, ]
    marker_censoring < -cbind(marker, other_markers)
    colnames(marker_censoring) < -paste("X", 1:ncol(marker_censoring), sep = "")
    fmla < - as.formula(paste("Surv(T,status)~", paste(paste("X", 1: ncol(
        marker_censoring), sep = ""), collapse = "+")))
    data_weight < -as.data.frame(cbind(data.frame(T=T, status=as.numeric(delta != 0)), marker_censoring))
    weights < - pec::ipcw(fmla, data=data_weight, method="cox", times=as.matrix(times), subjectTimes =
    data_weight[, "T"], subjectTimesLag = 1)
    }
    if (weighting == "aalen"){
    if (missing(other_markers)){marker_censoring < -marker}
    other_markers < -other_markers[order_T, ]
    marker_censoring < -cbind(marker, other_markers)
    colnames(marker_censoring) < -paste("X", 1:ncol(marker_censoring), sep = "")
    fmla < - as.formula(paste("Surv(T,status)~", paste(paste("X", 1: ncol(
        marker_censoring), sep = ""), collapse = "+")))
    data_weight < -as.data.frame(cbind(data.frame(T=T, status=as.numeric(delta != 0)), marker_censoring))
    weights < - pec::ipcw(fmla, data=data_weight, method="aalen", times=as.matrix(times), subjectTimes =
    data_weight[, "T"], subjectTimesLag = 1)
    }