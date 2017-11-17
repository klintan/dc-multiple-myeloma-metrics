##' ##-----------compare Kaplan-Meier to survival package---------
##'
##' dat <- SimSurv(30)
##' pfit <- prodlim(Surv(time,status)~1,data=dat)
##' pfit <- prodlim(Hist(time,status)~1,data=dat) ## same thing

def prodlim(formula, data, reverse=False):
    clustered = False
    if reverse:  ## estimation of censoring distribution
        model_type = 1
        if model_type == 1:
            if clustered:
                raise Exception('Not implemented yet')
            else:

                ## right censored not clustered
                # I think this one says: fit() function is corresponding to the call(); C is equivalent to call() similar to super() in python

                fit = {'name':'prodlimSRC', 'time': response["time"], 'status':response["status"], 0,
                    entrytime,caseweights, 0,N, 0, 0, NU,size.strata, 'time' = N, 'nrisk': N, 'nevent':N, 'ncens':
                    N, 'surv': N, 0, 'hazard': N, 'var.hazard': N, 'extra.double':
                    0, 'max.nc':0, 'ntimes': 1, 'ntimes.strata': NU, 'first.strata':
                    NU,reverse, 'model':0, 'independent':1, 'delayed':
                    delayed, 'weighted':weighted, 'PACKAGE':"prodlim"}
                # number of times
                NT = fit.ntimes

                Cout = {"time": fit.time[1:NT],
                        "n.risk" = fit.nrisk[1:NT],
                        "n.event" = fit.nevent[1:NT],
                       "n.lost" = fit.ncens[1:NT],
                        "surv" = fit.surv[1:NT],
                       "se.surv" = fit.surv[1:NT] * sqrt(pmax(0, fit$var.hazard[1:NT])),
                        "hazard" = fit.hazard[1:NT],
                        "first.strata" = fit.first.strata,
                         "size.strata" = fit.ntimes.strata,
                         "model" = "survival")
                Cout$maxtime < - max(Cout$time)
