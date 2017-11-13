import numpy as np
import functools
import pandas as pd
from ipcw import IPCW
from lifelines import KaplanMeierFitter
# import iid_decomposition as compute_iid_decomposition
import timeit


def reduce_concat(x, sep=""):
    return functools.reduce(lambda x, y: str(x) + sep + str(y), x)


def paste(*lists, sep=" ", collapse=None):
    result = map(lambda x: reduce_concat(x, sep=sep), zip(*lists))
    if collapse is not None:
        return reduce_concat(result, sep=collapse)
    return list(result)


def difftime():
    pass


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

    if len(times) < 1:
        raise ("Choose at least one time for computing the time-dependent AUC\n")

    if not weighting in ["marginal", "cox", "aalen"]:
        raise ("the weighting argument must be marginal (default), cox or aalen.\n")

    if weighting in ["cox", "aalen"] and isinstance(other_markers, (np.ndarray, np.generic)):
        raise ("argument other_markers must be numpy array/matrix\n")

    if weighting in ["cox", "aalen"] and other_markers:
        if np.shape(other_markers)[0] == len(marker):
            raise ("lengths of vector T, delta, marker and number of rows of other_markers have to be equal\n")

    # check if there are missing values, and delete rows with missing values
    # not implemented yet

    # Initialize variables
    start_computation_time = timeit.default_timer()

    n = len(T)
    n_marker = len(set(marker))  # get all unique markers ?
    n_times = len(times)
    if n_times == 1:
        times = [0] + times
        n_times = 2  # trick to use ipcw.cox() even if there is only one time

    times = sorted(times)

    # add possibility to use days or month in times
    times_names = ["t=" + str(x) for x in times]
    # times_names = paste(times, sep="t=")

    # Outputs
    AUC_1 = np.zeros(n_times)  # I don't think these are necessary but will stay for now
    AUC_2 = np.zeros(n_times)  # I don't think these are necessary but will stay for now
    CumInci = np.zeros(n_times)
    surv = np.zeros(n_times)

    # names(AUC_1) < -times_names
    # names(AUC_2) < -times_names
    # names(CumInci) < -times_names
    # names(surv) < -times_names


    stats_data = np.zeros((n_times, 4))
    Stats = pd.DataFrame(stats_data, columns=["Cases", "survivor at t", "Other events at t", "Censored at t"],
                         index=times_names)
    # rownames(Stats) < -times_names # set the index to the times names?

    # computation of weights (1/2)

    # we need to order to use the pec::ipcw() function
    indices = np.argsort(T)
    T = np.array(sorted(T))
    delta = delta.values[indices]
    marker = marker.values[indices]
    weights = {}
    # use ipcw function from pec package
    if (weighting == "marginal"):
        Surv = "formula"  # Surv(failure_time, status)
        ipcw_data = np.array((T, delta)).transpose()
        df_ipcw_data = pd.DataFrame(data=ipcw_data, index=ipcw_data[:, 0], columns=['failure_time', 'status'])
        ipcw = IPCW(formula=Surv, data=df_ipcw_data, method="marginal", times=times, subjectTimes=T, what=["IPCW.times", "IPCW.subject.times"], subjectTimesLag=1)
        weights = ipcw.fit()

    if (weighting == "cox"):
        raise ("Cox weighing not yet supported")
    if (weighting == "aalen"):
        raise ("Aalen weighing not yet supported")

    # we order by marker values (in order to compute Se and Sp)
    order_marker = sorted(marker)
    # concat by column all of them using order marker (to get corresponding values in order)
    Mat_data = pd.DataFrame([T, delta, marker][order_marker,], columns=["T", "delta",
                                                                        "marker"])  # all matrices should be defined as Pandas dataframes instead

    # Create some weights
    Weights_cases_all = 1 / (weights.subjectTimes * n)
    Weights_cases_all = Weights_cases_all[order_marker]

    #  Make TP and FP outputs if needed
    if ROC == True:
        FP_1 = np.matrix(np.zeros(n_marker + 1), np.zeros(n_times))
        TP = np.matrix(np.zeros(n_marker + 1), np.zeros(n_times))
        FP_2 = np.matrix(np.zeros(n_marker + 1), np.zeros(n_times))
    else:
        FP_1 = None
        FP_2 = None
        TP = None

    # loop on all timepoints t
    for t in range(1, n_times):
        Cases = Mat_data["T"] < times[t] & Mat_data["delta"] == cause
        Controls_1 = Mat_data["T"] > times[t]
        Controls_2 = Mat_data["T"] < times[t] & Mat_data["delta"] != cause & Mat_data["delta"] != 0
        if weighting != "marginal":
            Weights_controls_1 = 1 / (weights.times(t) * n)
        else:
            # Weights_controls_1 = 1 / (weights.times[t] * n), times=n
            pass

        Weights_controls_1 = Weights_controls_1[order_marker]

        Weights_cases = Weights_cases_all
        Weights_controls_2 = Weights_cases_all

        Weights_cases[Cases] = 0  # DF is.na ?
        Weights_controls_1[Controls_1] = 0  # DF is.na ?
        Weights_controls_2[Controls_2] = 0  # DF is.na ?

        den_TP_t = sum(Weights_cases)
        den_FP_1_t = sum(Weights_controls_1)
        den_FP_2_t = sum(Weights_controls_2) + sum(Weights_controls_1)

        if den_TP_t != 0:
            TP_tbis = np.array(0, np.cumsum(Weights_cases)) / den_TP_t
            TP_t = TP_tbis[pd.duplicated(marker[order_marker])]  # get all not duplicated (needs to be a dataframe)
        else:
            TP_t = None

        if (den_FP_1_t != 0):
            FP_1_tbis = np.array(0, np.cumsum(Weights_controls_1)) / den_FP_1_t
            FP_1_t = FP_1_tbis[pd.duplicated(marker[order_marker])]
        else:
            FP_1_t = None

        if den_FP_2_t != 0:
            FP_2_tbis = np.array(0, np.cumsum(Weights_controls_1) + np.cumsum(Weights_controls_2)) / den_FP_2_t
            FP_2_t = FP_2_tbis[pd.duplicated(marker[order_marker])]
        else:
            FP_2_t = None

        # internal function to compute an area under a curve by trapezoidal rule
        def airetrap(Abs, Ord):
            nobs = len(Abs)
            dAbs = Abs[-1] - Abs[-nobs]
            mil = (Ord[-nobs] + Ord[-1]) / 2
            area = sum(dAbs * mil)
            return (area)

        if den_TP_t * den_FP_1_t != 0:
            AUC_1[t] = airetrap(FP_1_t, TP_t)
        else:
            AUC_1[t] = None

        if (den_TP_t * den_FP_2_t != 0):
            AUC_2[t] = airetrap(FP_2_t, TP_t)
        else:
            AUC_2[t] = None

        if ROC == True:
            TP[:, t] = TP_t
            FP_1[:, t] = FP_1_t
            FP_2[:, t] = FP_2_t

    CumInci[t] = np.array(den_TP_t)
    surv[t] = np.array(den_FP_1_t)
    Stats[t, :] = np.array(sum(Cases), sum(Controls_1), sum(Controls_2),
                           n - sum(Cases) - sum(Controls_1) - sum(Controls_2))

    inference = None
    if iid == True:
        if weighting != "marginal":
            raise (
                "Error : Weighting must be marginal for computing the iid representation \n Choose iid=FALSE or weighting=marginal in the input arguments")

    else:
        # create iid representation required for inference procedures
        out_iid = np.array("list", n_times)
        names[out_iid] = paste("t=", times, sep="")
        vect_iid_comp_time = rep(NA, times=n_times)
        names[vect_iid_comp_time] = paste("t=", times, sep="")
        mat_iid_rep = np.matrix(np.zeros(n), np.zeros(n_times))
        colnames[mat_iid_rep] = paste("t=", times, sep="")
        mat_iid_rep_star = np.matrix(np.zeros(n), np.zeros(n_times))
        colnames[mat_iid_rep_star] = paste("t=", times, sep="")
        vetc_se = rep(NA, times=n_times)
        names[vetc_se] = paste("t=", times, sep="")
        vetc_sestar = rep(NA, times=n_times)
        names[vetc_sestar] = paste("t=", times, sep="")

    # compute iid for Kaplan Meier
    kmf = KaplanMeierFitter()
    MatInt0TcidhatMCksurEff = kmf.fit(T, event_observed=delta)
    # MatInt0TcidhatMCksurEff = Compute.iid.KM(times=T, status=delta)
    for j in range(1, n_times):
        # compute iid representation when AUC can be computed
        if AUC_1[j] or AUC_2[j]:
            out_iid[[j]] = compute_iid_decomposition(weights, T, delta, marker, t=times[j], n=n, cause=cause,
                                                     F01t=CumInci[j], St=surv[j],
                                                     MatInt0TcidhatMCksurEff=MatInt0TcidhatMCksurEff)
        else:
            out_iid[[j]] = None
        # browser()
        # save output for inference for AUC_1 when AUC_1 can be computed
        if AUC_1[j]:
            mat_iid_rep_star[:, j] = out_iid[j].iid_representation_AUCstar
        vetc_sestar[j] = out_iid[j].seAUCstar
        vect_iid_comp_time[j] = out_iid[j].computation_times

        # save output for inference for AUC_2 when AUC_2 can be computed
        if AUC_2[j]:
            mat_iid_rep[:, j] = out_iid[j].iid_representation_AUC
        vetc_se[j] = out_iid[j].seAUC
        vect_iid_comp_time[j] = out_iid[[j]].computation_times

        inference = {'mat_iid_rep_2': mat_iid_rep,
                     'mat_iid_rep_1': mat_iid_rep_star,
                     'vect_sd_1': vetc_sestar,
                     'vect_sd_2': vetc_se,
                     'vect_iid_comp_time': vect_iid_comp_time
                     }

    # output if there is competing risks or not
    if max(Stats[:, 3]) == 0:
        out = {'TP': TP, 'FP': FP_1, 'AUC': AUC_1, 'times': times, 'CumulativeIncidence': CumInci, 'survProb': surv,
               'n': n, 'Stats': Stats[:, np.array(1, 2, 4)], 'weights': weights,
               'inference': inference,
               'computation_time': difftime(stop_computation_time, start_computation_time, units="secs"), 'iid': iid}
        # class(out) < - "ipcwsurvivalROC"
        print(out)
    else:
        out = {'TP': TP, 'FP_1': FP_1, 'AUC_1': AUC_1, 'FP_2': FP_2, 'AUC_2': AUC_2, 'times': times,
               'CumulativeIncidence': CumInci, 'survProb': surv, 'n': n, 'Stats': Stats, 'weights': weights,
               'inference': inference,
               'computation_time': difftime(stop_computation_time, start_computation_time, units="secs"), 'iid': iid}
        # class(out) < - "ipcwcompetingrisksROC"
        print(out)
