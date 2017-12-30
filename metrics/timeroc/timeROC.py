import numpy as np
import pandas as pd
from metrics.ipcw.ipcw import IPCW
from lifelines import KaplanMeierFitter
# import iid_decomposition as compute_iid_decomposition
import timeit
from sklearn.metrics import auc


def calculate_weights_cases_all(mat_data, order, n):
    # weights for all cases
    weights_cases_all = 1 / (mat_data['T'].values * n)
    # sort using arg_sort list
    return weights_cases_all[order]


def case_masking(mat_data, times, t, cause):
    times_t = mat_data["T"] < times[t]
    delta_cause = mat_data["delta"] == cause
    delta_not_cause = mat_data["delta"] != cause
    delta_one = mat_data["delta"] != 0

    cases = mat_data[times_t & delta_cause]
    controls_1 = mat_data[times_t]
    controls_2 = mat_data[times_t & delta_not_cause & delta_one]

    return (cases, controls_1, controls_2)

def matrix_transpose(*args):
    return np.array((args)).transpose()


def timeROC(T, delta, marker, cause, times, other_markers=None, weighting="marginal", ROC=True, iid=False):
    '''
    A function which returns time dependent ROC values.
    Adapted from R timeROC https://github.com/cran/timeROC/blob/master/R/timeROC_3.R

    :param T: vector of observed failure times
    :param delta: vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...)
    :param marker: vector of marker values
    :param cause: the value that indicates the main event of interest
    :param times: vector of times you want to compute the time dependent AUC.
    :param other_markers: (default is NULL, should be a matrix) other markers that can be associated with the censoring mechanism
    :param weighting: (default is "marginal") weighting technique for IPCW : "marginal" for Kaplan-Meier, "cox" for proportional hazards cox model, "aalen" for additive aalen model
    :param ROC: if TRUE, then save True Positive fraction (Sensitivity) and False Positive fraction (1-Specificity) for all time in vetor times
    :param iid: TRUE or FALSE, indicates if we want to compute the iid representation of the AUC estimator
    :return:
        TP : if ROC=TRUE, the matrix of True Positive fraction (sensitivty),
        with  columns correspond to time of vector (ordered) times
        FP_1 : if ROC=TRUE, the matrix of False Positive fraction (1-Specificity),
        where a "control" is defined as an event-free subject at time t,
        with  columns correspond to time of vector (ordered) times
        FP_2 : if ROC=TRUE, the matrix of False Positive fraction (1-Specificity),
        where a "control" is defined as any subject that is not a case,
        with  columns correspond to time of vector (ordered) times
        AUC_1 : vector of AUC for each times in vector (ordered) times computed with FP_1
        AUC_2 : vector of AUC for each times in vector (ordered) times computed with FP_2
        times : the input vector (ordered) times. If there is only one time in the input,
        then 0 is added
        survProb: the vector of probabilities to be event free (all causes) at each time in
        vector (ordered) times (denominator of False positive fraction)
        CumulativeIncidence : the vector of the cumulative incidences of the type of event
        of interest at each time in vector (ordered) times
        n : number of observations
        Stats : matrix with numbers of observed cases and controls (both definitions), and
        censored before time point, for each  time in vector (ordered) times
        weights : an object of class "ipcw" (see "pec" package) with all information about the weights
        computation_time : the computation time

    '''

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
    n_marker = len(set(marker))  # number of unique markers
    n_times = len(times)

    if n_times == 1:
        times = [0] + times
        n_times = 2  # trick to use ipcw.cox() even if there is only one time

    times = sorted(times)

    # add possibility to use days or month in times
    times_names = ["t=" + str(x) for x in times]

    # Outputs
    AUC_1 = np.zeros(n_times)  # I don't think these are necessary but will stay for now
    AUC_2 = np.zeros(n_times)  # I don't think these are necessary but will stay for now
    CumInci = np.zeros(n_times)
    surv = np.zeros(n_times)

    stats_data = np.zeros((n_times, 4))
    Stats = pd.DataFrame(stats_data, columns=["Cases", "survivor at t", "Other events at t", "Censored at t"],
                         index=times_names)

    indices = np.argsort(T)

    T = np.array(sorted(T))
    delta = delta.values[indices]
    marker = marker.values[indices]
    weights = {}
    IPCW_data = matrix_transpose(T, delta)
    temp_data = matrix_transpose(T, delta, marker)

    # use ipcw function from pec package
    if weighting == "marginal":
        df_ipcw_data = pd.DataFrame(data=IPCW_data, index=IPCW_data[:, 0], columns=['failure_time', 'status'])
        ipcw = IPCW(formula=None, data=df_ipcw_data, method="marginal", times=times, subjectTimes=T,
                    what=["IPCW.times", "IPCW.subject.times"], subjectTimesLag=1)
        # TODO this should be ipcw.fit()
        weights = ipcw.marginal()

    if weighting == "cox":
        raise NotImplementedError('Cox weighing not yet supported')

    if weighting == "aalen":
        raise NotImplementedError('Aalen weighing not yet supported')

    # we order by marker values (in order to compute Se and Sp)
    order_marker = np.argsort(marker)

    # concat by column all of them using order marker (to get corresponding values in order)
    Mat_data = pd.DataFrame(temp_data, columns=["T", "delta", "marker"])
    Mat_data.sort_values("marker", ascending=False, inplace=True)

    Weights_cases_all = calculate_weights_cases_all(Mat_data, order_marker, n)

    #  Make TP and FP outputs if needed
    if ROC:
        FP_1 = np.zeros((n_marker, n_times))
        TP = np.zeros((n_marker, n_times))
        FP_2 = np.zeros((n_marker, n_times))
    else:
        FP_1 = None
        FP_2 = None
        TP = None

    # loop on all timepoints t
    cause = float(cause)
    for t in range(0, n_times):
        Cases, Controls_1, Controls_2 = case_masking(Mat_data, times, t, cause)

        if weighting != "marginal":
            Weights_controls_1 = 1 / (times[t] * n)
        else:
            Weights_controls_1 = np.repeat(1 / (times[t] * n), n)

        # Sort, skip for now, probably unnecessary
        Weights_controls_1 = Weights_controls_1[order_marker]

        Weights_cases = Weights_cases_all
        Weights_controls_2 = Weights_cases_all
        # All the patients in Cases should be 0 for Weights_cases
        Weights_cases[Cases.index] = 0
        # All the patients in Weights_controls_1 should be 0 for Controls_1
        Weights_controls_1[Controls_1.index] = 0  # DF is.na ?
        # All the patients in Weights_controls_2 should be 0 in Controls_2 ?
        Weights_controls_2[Controls_2.index] = 0  # DF is.na ?

        den_TP_t = sum(Weights_cases)
        den_FP_1_t = sum(Weights_controls_1)
        den_FP_2_t = sum(Weights_controls_2) + sum(Weights_controls_1)

        if den_TP_t != 0:
            TP_tbis = np.cumsum(Weights_cases) / den_TP_t
            TP_t = TP_tbis[np.unique(marker[order_marker], return_index=True)[1]]
        else:
            TP_t = None

        if den_FP_1_t != 0:
            FP_1_tbis = np.cumsum(Weights_controls_1) / den_FP_1_t
            FP_1_t = FP_1_tbis[np.unique(marker[order_marker], return_index=True)[1]]
        else:
            FP_1_t = None

        if den_FP_2_t != 0:
            FP_2_tbis = (np.cumsum(Weights_controls_1) + np.cumsum(Weights_controls_2)) / den_FP_2_t
            FP_2_t = FP_2_tbis[np.unique(marker[order_marker], return_index=True)[1]]
        else:
            FP_2_t = None

        if den_TP_t * den_FP_1_t != 0:
            # sorting that sklearn auc uses
            # order = np.lexsort((FP_1_t, TP_t))
            # FP_1_t, TP_t = FP_1_t[order], TP_t[order]
            # it needs to be ordered, it should probably be ordered already
            AUC_1[t] = auc(FP_1_t, TP_t)

        else:
            AUC_1[t] = None

        if den_TP_t * den_FP_2_t != 0:
            # it needs to be ordered, it should probably be ordered already
            AUC_2[t] = auc(FP_2_t, TP_t)

        else:
            AUC_2[t] = None

        if ROC == True:
            TP[:, t] = TP_t
            FP_1[:, t] = FP_1_t
            FP_2[:, t] = FP_2_t

        CumInci[t] = den_TP_t
        surv[t] = den_FP_1_t
        # Stats.iloc[t] = [sum(Cases.values), sum(Controls_1.values), sum(Controls_2.values), n - sum(Cases.values) - sum(Controls_1.values) - sum(Controls_2.values)]

        # Stats[t, :] = np.array([sum(Cases), sum(Controls_1), sum(Controls_2),
        #                       n - sum(Cases) - sum(Controls_1) - sum(Controls_2)])

    inference = None
    if iid == True:
        if weighting != "marginal":
            raise (
                "Error : Weighting must be marginal for computing the iid representation \n Choose iid=FALSE or weighting=marginal in the input arguments")

        else:
            # create iid representation required for inference procedures
            out_iid = np.zeros(n_times)
            # names[out_iid] = paste("t=", times, sep="")

            vect_iid_comp_time = np.zeros(n_times)
            # names[vect_iid_comp_time] = paste("t=", times, sep="")

            mat_iid_rep = np.zeros(n, n_times)
            # colnames[mat_iid_rep] = paste("t=", times, sep="")

            mat_iid_rep_star = np.zeros(n, n_times)
            # colnames[mat_iid_rep_star] = paste("t=", times, sep="")

            vetc_se = np.zeros(n_times)
            # names[vetc_se] = paste("t=", times, sep="")

            vetc_sestar = np.zeros(n_times)
            # names[vetc_sestar] = paste("t=", times, sep="")

        # compute iid for Kaplan Meier
        kmf = KaplanMeierFitter()
        MatInt0TcidhatMCksurEff = kmf.fit(T, event_observed=delta)
        # MatInt0TcidhatMCksurEff = Compute.iid.KM(times=T, status=delta)
        for j in range(1, n_times):
            # compute iid representation when AUC can be computed
            if AUC_1[j] or AUC_2[j]:
                raise NotImplementedError("IID decomposition not implemented yet")
                # out_iid[j] = compute_iid_decomposition(weights, T, delta, marker, t=times[j], n=n, cause=cause,
                #                                         F01t=CumInci[j], St=surv[j],
                #                                         MatInt0TcidhatMCksurEff=MatInt0TcidhatMCksurEff)
            else:
                out_iid[j] = None
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

    stop_computation_time = timeit.default_timer()

    # output if there is competing risks or not
    if max(Stats['Censored at t'].values) == 0:
        out = {'TP': TP, 'FP': FP_1, 'AUC': AUC_1, 'times': times, 'CumulativeIncidence': CumInci, 'survProb': surv,
               'n': n, 'Stats': Stats[['Cases', 'Other events at t', 'Censored at t']], 'weights': weights,
               'inference': inference,
               'computation_time': stop_computation_time - start_computation_time, 'iid': iid}
        # class(out) < - "ipcwsurvivalROC"
    else:
        out = {'TP': TP, 'FP_1': FP_1, 'AUC_1': AUC_1, 'FP_2': FP_2, 'AUC_2': AUC_2, 'times': times,
               'CumulativeIncidence': CumInci, 'survProb': surv, 'n': n, 'Stats': Stats, 'weights': weights,
               'inference': inference,
               'computation_time': stop_computation_time - start_computation_time, 'iid': iid}
        # class(out) < - "ipcwcompetingrisksROC"
    return out
