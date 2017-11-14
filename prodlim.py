##' ##-----------compare Kaplan-Meier to survival package---------
##'
##' dat <- SimSurv(30)
##' pfit <- prodlim(Surv(time,status)~1,data=dat)
##' pfit <- prodlim(Hist(time,status)~1,data=dat) ## same thing

def prodlim(formula,data, reverse=False):
    clustered = False
    if reverse:  ## estimation of censoring distribution
        model_type = 1
        if model_type == 1:
            if clustered:
                raise Exception('Not implemented yet')
            else:

                ## right censored not clustered
                # I think this one says: fit() function is corresponding to the call(); C is equivalent to call() similar to super() in python

                fit  = .C("prodlimSRC",as.double(response[, "time"]), as.double(response[, "status"]), integer(0),as.double(
                    entrytime),as.double(caseweights), integer(0),as.integer(N), integer(0), integer(0),as.integer(
                    NU),as.integer(size.strata), time = double(N), nrisk = double(N), nevent = double(N), ncens = double(
                    N), surv = double(N), double(0), hazard = double(N), var.hazard = double(N), extra.double = double(
                    0), max.nc = integer(0), ntimes = integer(1), ntimes.strata = integer(NU), first.strata = integer(
                    NU),as.integer(reverse), model =as.integer(0), independent =as.integer(1), delayed =as.integer(
                    delayed), weighted =as.integer(weighted), PACKAGE = "prodlim")
                # number of times
                NT =  fit$ntimes

                Cout = {"time": fit$time[1:NT],
                                           "n.risk" = fit$nrisk[1:NT],
                                                          "n.event" = fit$nevent[1:NT],
                                                                          "n.lost" = fit$ncens[1:NT],
                                                                                         "surv" = fit$surv[1:NT],
                                                                                                      "se.surv" = fit$surv[
                                                                                                                      1:NT] * sqrt(
                    pmax(0, fit$var.hazard[1:NT])),
                "hazard" = fit$hazard[1:NT],
                               "first.strata" = fit$first.strata,
                                                    "size.strata" = fit$ntimes.strata,
                                                                        "model" = "survival")
                Cout$maxtime < - max(Cout$time)
