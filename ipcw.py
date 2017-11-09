# Inverse of the probability of censoring weights (IPCW)
# usually refer to the probabilities of not being censored at
# certain time points. These probabilities are also the
# values of the conditional survival function of the censoring
# time given covariates. The function ipcw estimates the conditional
# survival function of the censoring times and derives the weights.
#IMPORTANT: the data set should be ordered, order(time,-status) in order
# to get the values IPCW.subjectTimes in the right order for some choices of method.

import pandas as pd
import lifelines

def ipcw(formula, data, method, args, times, subjectTimes, what, subjectTimesLag = 1,):
    pass