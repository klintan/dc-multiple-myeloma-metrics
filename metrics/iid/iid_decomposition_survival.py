import timeit



def difftime():
    pass


def compute_iid_decomposition_survival(t, n, cause, F01t, St, weights, T, delta, marker, MatInt0TcidhatMCksurEff):
    """
    # CAUTION : T,delta,marker,weights should be order by order(T)

    :param t: time at which we compute the i.i.d decomposition (Influence function)
    :param n: sample size
    :param cause: the value that indicates the main event of interest
    :param F01t: the cumulative incidence function of the main event at time t
    :param St: object weights, output of the main function, order by order(T) (by default with ipcw() function of package pec)
    :param weights: vector of observed failure times, order by order(T)
    :param T: vector of indicator of status (0 for censoring, 1 for type of event one, 2 for type of event two and so on...),order by order(T)
    :param delta: vector ofmarker values,order by order(T)
    :param marker: vector of times for wich we compute the AUCs
    :param MatInt0TcidhatMCksurEff:
    :return:
    """

    start_total = timeit.default_timer()

    # indicator vectors
    Cases = (T < t & delta == cause)
    Controls_1 = (T > t)
    # vectors which indicates the indexes of Cases and the Controls
    which_Cases = which(T < t & delta == cause)
    which_Controls_1 = which(T > t)
    # compute the weights.
    Weights_cases_all = 1 / (weights$IPCW.subjectTimes * n)
    Weights_cases = Weights_cases_all
    Weights_cases[Cases] = 0  # (0 if not a case)
    Weights_controls_1 = rep(1 / (weights.IPCW.times[which(weights.times == t)] * n), times=n)
    Weights_controls_1[Controls_1] = 0  # (0 if not a control)
    # compute vector indicator of censoring (event is censoring !)
    indic_Cens =as.numeric(delta == 0)
    # compute the matrix with all information. The matrix is order by order(t)
    Mat_data = cbind(T, delta, indic_Cens, marker, Cases, Controls_1, Weights_cases, Weights_controls_1)

    ## MatInt0TcidhatMCksurEff = Compute.iid.KM(times=T,status=delta)
    # {{{ STEP : Compute terms  {\hat{h}_{tij}}_1 and {\hat{h}_{tij}}_2
    # start_htij=Sys.time()
    # function that eats the matrix W1 (defined just after) that depends on subject i and returns
    # the vector of {\hat{h}_{tij}}_1
    def htij1(V, tps=t):
        pass
        #return as.numeric(V[, 1] > tps)*(as.numeric(V[, 4] > V[, 2]) + 0.5 *as.numeric(V[, 4] == V[, 2])) *(V[, 3]*V[, 5])*(n * n)

    # compute frequencies of cases and controls to define
    # the size of the matrix  Mathtij1
    nb_Cases = sum(T < t & delta == cause)
    nb_Controls_1 = sum(T > t)
    # To save computation time, we loop only on control 1 for Mathtij1
    Mat_data_cont1 = Mat_data[which_Controls_1,]
    # initialise  Mathtij1  with its right size !
    Mathtij1 = matrix(NA, nb_Controls_1, nb_Cases)
    # loop on all cases i. We loop only on Cases to save computation time !
    for i in range(which_Cases):
        #W1 = cbind(Mat_data_cont1[, c("T", "marker")], rep(Mat_data[i, c("Weights_cases")], nb_Controls_1), rep(Mat_data[i, c("marker")], nb_Controls_1), Mat_data_cont1[, c("Weights_controls_1")])
        W1 = {}
        # fill the column i of  Mathtij1 and  Mathtij2
        #Mathtij1[, which(i == which_Cases)] = htij1(W1)
        Mathtij1 = {}
    # matrix Mathtij1  : i for columns, j for rows
    # browser() # nice function for debugging !
    # stop_htij=Sys.time()
    # print(difftime(stop_htij,start_htij,units="sec"))
    # compute \hat{h}_t
    ht = (sum(Mathtij1)) / (n * n)
    # vector of \hat{f}_{i1t}
    #vect_dit =as.numeric(Mat_data[, c("T")] <= t) * as.numeric(Mat_data[, c("delta")] == cause)*Mat_data[, c("Weights_cases")]*n
    vect_dit = {}
    # We can check we have F01t by mean(vect_dit)
    # print("F01t ??")
    # print(c(mean(vect_dit),F01t))

    #  Final step : to compute iid representation of AUC^*(t)
    start_iid_AUC1 = Stimeit.default_timer()
    # Compute the vecor of all sum_{i=1}^n of {\hat{h}_{tij}}_1 for all j
    #colSums_Mathtij1 = rep(0, n)  # initialise at 0
    colSums_Mathtij1 = {}
    colSums_Mathtij1[which_Cases] = colSums(Mathtij1)  # when i is a case,  then we sum the column of  Mathtij1
    # Compute the vecor of all sum_{j=1}^n of {\hat{h}_{tij}}_1 for all i
    rowSums_Mathtij1 = rep(0, n)  # initialize at 0
    rowSums_Mathtij1[which_Controls_1] = rowSums(Mathtij1)  # when  j is a control 1, then we sum the row of  Mathtij1
    hathtstar = (sum(Mathtij1)) / (n * n)
    # print("AUC1 ???")
    # print(hathtstar/(F01t*St))
    # compute the vector of \frac{1_{\tilde{T}_i>=t}}{ \hat{S}_{\tilde{T}}(t)}
    vect_Tisupt = as.numeric(Mat_data[, c("T")] > t) / (sum(as.numeric(Mat_data[, c("T")] > t)) / n )
    def sum_ij_a_k_fixe(k):
        Pour_sum_ij_a_k_fixe = t(Mathtij1) * (1 + MatInt0TcidhatMCksurEff[which_Cases, k])
        Pour_sum_ij_a_k_fixe_3 = vect_dit * (1 + MatInt0TcidhatMCksurEff[, k])
        Pour_sum_ij_a_k_fixe_3b = (hathtstar) * (vect_Tisupt + (1 / F01t) * (Pour_sum_ij_a_k_fixe_3 - F01t))
        La_sum_ij_a_k_fixe = sum(Pour_sum_ij_a_k_fixe) / n - sum(Pour_sum_ij_a_k_fixe_3b)
        return (La_sum_ij_a_k_fixe)

    # print("F01t*St")
    # print(F01t*St)
    Les_sum_ij_a_k_fixe = (sapply(1:n, sum_ij_a_k_fixe)) / (F01t * St)
    Les_sum_ik_a_j_fixe = (rowSums_Mathtij1 - n * hathtstar) / (F01t * St)
    Les_sum_jk_a_i_fixe = (colSums_Mathtij1 - n * hathtstar * (vect_Tisupt + (1 / F01t) * (vect_dit - F01t))) / (
        F01t * St)
    # We compute the iid representation of the AUC estimator
    hatIFstar = (Les_sum_ij_a_k_fixe + Les_sum_ik_a_j_fixe + Les_sum_jk_a_i_fixe) / (n)
    stop_iid_AUC1 = Sys.time()
    # }}}
    # we compute the standard error of the AUC estimator
    seAUCstar = sd(hatIFstar) / sqrt(n)
    # browser() # nice function for debugging
    stop_total = timeit.default_timer()
    total_time = difftime(stop_total, start_total, units="secs")
    total_time_iid_AUC1 = difftime(stop_iid_AUC1, start_iid_AUC1, units="secs")
    computation_times = c(total_time)
    names[computation_times] = np.array("total_time")
    return {'iid_representation_AUC': rep(NA, n),
            'iid_representation_AUCstar': hatIFstar,
            'seAUC': None, 'seAUCstar': seAUCstar,
            'computation_times': computation_times}
