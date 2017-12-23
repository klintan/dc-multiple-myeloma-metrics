import numpy as np
from itertools import compress

def IntegrateAUC(AUC, utimes, St, tmax, weight = "rescale"):
    '''
    A function which integrates the AUC over times utimes
    :param AUC: AUC values for each times
    :param utimes: the time values specified
    :param St: tempAUC$survProb
    :param tmax: the max(times)
    :param weight:
    :return:
    '''
    if weight not in ['conditional', 'rescale']:
        raise Exception("error in weight choice")

    ft = np.zeros(len(St))
    ft[1] = 1 - St[1]

    for j in range(2,len(St)):
        ft[j] = St[j - 1] - St[j]

    #filter list using bool list
    mIndex = len(list(compress(utimes, utimes <= tmax)))
    www = 2 * ft * St
    wTotal = sum(www[1:mIndex])
    # Smax <- St[min(mIndex + 1, length(St))]
    Smax = St[min([mIndex, len(St)-1])]
    if weight == "conditional":
        w = 2 * ft * St/wTotal
    else:
        w = 2 * ft * (St - Smax)/((1 - Smax)**2)

    iAUC = sum(w[1:mIndex] * AUC[1:mIndex])
    return(iAUC)
