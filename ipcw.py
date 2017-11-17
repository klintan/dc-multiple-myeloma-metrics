import pandas as pd
from lifelines import KaplanMeierFitter
import numpy as np
import statsmodels.formula.api as smf


class IPCW():
    """
    Estimation of censoring probabilities

    This function is used internally to obtain
    inverse of the probability of censoring weights.

    Inverse of the probability of censoring weights (IPCW) usually refer to the
    probabilities of not being censored at certain time points. These
    probabilities are also the values of the conditional survival function of
    the censoring time given covariates. The function ipcw estimates the
    conditional survival function of the censoring times and derives the
    weights.

    IMPORTANT: the data set should be ordered, \code{order(time,-status)} in
    order to get the values \code{IPCW.subject.times} in the right order for some
    choices of \code{method}.

    @aliases ipcw ipcw.none ipcw.marginal ipcw.nonpar ipcw.cox ipcw.aalen
    @param formula A survival formula like, \code{Surv(time,status)~1}, where
    as usual status=0 means censored. The status variable is internally
      reversed for estimation of censoring rather than survival
      probabilities. Some of the available models (see argument
     \code{model}) will use predictors on the right hand side of the
     formula.
     @param data The data used for fitting the censoring model
     @param method Censoring model used for estimation of the
     (conditional) censoring distribution.
     @param args A list of arguments which is passed to method
     @param times For \code{what="IPCW.times"} a vector of times at
     which to compute the probabilities of not being censored.
     @param subject.times For \code{what="IPCW.subject.times"} a vector of
     individual times at which the probabilities of not being censored
     are computed.
     @param lag If equal to \code{1} then obtain
     \code{G(T_i-|X_i)}, if equal to \code{0} estimate the conditional
     censoring distribution at the subject.times,
     i.e. (\code{G(T_i|X_i)}).
     @param what Decide about what to do: If equal to
     \code{"IPCW.times"} then weights are estimated at given
     \code{times}.  If equal to \code{"IPCW.subject.times"} then weights
     are estimated at individual \code{subject.times}.  If missing then
     produce both.
     @param keep Which elements to add to the output. Any subset of the vector \code{c("times","fit","call")}.
     @return A list with elements depending on argument \code{keep}. \item{times}{The times at which weights are estimated}
     \item{IPCW.times}{Estimated weights at \code{times}}
     \item{IPCW.subject.times}{Estimated weights at individual time values
     \code{subject.times}} \item{fit}{The fitted censoring model}
     \item{method}{The method for modelling the censoring distribution}
     \item{call}{The call}
     @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
     @keywords survival
     @examples

     library(prodlim)
     library(rms)
     dat=SimSurv(30)

     dat <- dat[order(dat$time),]

     # using the marginal Kaplan-Meier for the censoring times

     WKM=ipcw(Hist(time,status)~X2,
       data=dat,
       method="marginal",
       times=sort(unique(dat$time)),
       subject.times=dat$time,keep=c("fit"))
     plot(WKM$fit)
     WKM$fit

     # using the Cox model for the censoring times given X2
     library(survival)
     WCox=ipcw(Hist(time=time,event=status)~X2,
       data=dat,
       method="cox",
       times=sort(unique(dat$time)),
       subject.times=dat$time,keep=c("fit"))
     WCox$fit

     plot(WKM$fit)
     lines(sort(unique(dat$time)),
           1-WCox$IPCW.times[1,],
           type="l",
           col=2,
           lty=3,
           lwd=3)
     lines(sort(unique(dat$time)),
           1-WCox$IPCW.times[5,],
           type="l",
           col=3,
           lty=3,
           lwd=3)

     # using the stratified Kaplan-Meier
     # for the censoring times given X2

     WKM2=ipcw(Hist(time,status)~X2,
       data=dat,
       method="nonpar",
       times=sort(unique(dat$time)),
       subject.times=dat$time,keep=c("fit"))
     plot(WKM2$fit,add=FALSE)
    """

    def __init__(self, formula, data, method, times, subjectTimes, what=None, args=None, subjectTimesLag=1):
        if not what:
            raise Exception('Missing what: will produce IPCW.times and IPCW.subject.times.')
        #if not what or not what in ["IPCW.times", "IPCW.subject.times"]:
        #    raise Exception(len(times) > 0)

        self.formula = formula
        self.data = data
        self.method = method
        self.args = args
        self.times = times
        self.subject_times = subjectTimes
        self.what = what
        self.lag = subjectTimesLag
        self.keep = None
        self.lag = None
        self.methods = {'none': self.none(), 'rfsrc': self.rfsrc(), 'forest': self.forest(), 'marginal': self.marginal(),
                        'nonpar': self.nonpar(), 'cox': self.cox()}

        # dummy variables for now
        self.call = {}
        self.fit = {}
        self.kmf = KaplanMeierFitter()

        # check input arguments
        self.args_check()

    def fit(self):
        run_method = self.methods[self.method]
        run_method(self)

    @staticmethod
    def output(out, keep, times, fit, call):
        if "call" in keep:
            output = c(out, list(call=call))
        if "times" in keep:
            output = c(out, list(times=times))
        if "fit" in keep:
            output = c(out, list(fit=fit))
        return output

    def args_check(self):
        if not self.lag:
            self.lag = 1
        if not self.what:
            self.what = ["IPCW.times", "IPCW.subject.times"]

    def none(self):
        """
        Uncensored data: return just array of 1:s
        :return:
        """

        #  weights at requested times
        if "IPCW.times" in self.what:
            self.times = np.ones(len(self.times))
        else:
            self.times = None

        # weights at subject specific event times
        if "IPCW.subject.times" in self.what:
            self.subject_times = np.ones(len(self.subject_times))
        else:
            self.subject_times = None

        return {'times': self.times, 'subject_times': self.subject_times, 'method': self.method, 'call': None}

    # reverse Random Survival Forests
    def rfsrc(self):
        # call = match.call()  ## needed for refit in crossvalidation loop
        EHF = {}
        # EHF = prodlim.EventHistory.frame(formula, data,specials=None, unspecialsDesign=False)
        # wdata = data.frame(cbind(unclass(EHF.event.history), EHF.design))
        # wdata <- as.data.frame(EHF)
        wdata = pd.DataFrame(EHF)
        #wdata.status = 1 - wdata.status

        # R function update: e update function to start with a formula from a model that we have already fitted and to
        #  specify the terms that we want to add or remove (we need to use python statsmodel)
        # wform = update(formula, "Surv(time,status)~.")
        ## require(randomForestSRC)
        # if not (NROW(na.omit(wdata)) > 0):
        #    raise("stop")
        if not self.args:
            self.args = {'ntree':1000, 'importance': None, 'bootstrap':None, 'nodesize': len(self.data.index)/2}

            # fit = do.call(randomForestSRC.rfsrc, c(list(wform, data=wdata), args))
            # fit.call = None
            #  weights at requested times
            #  predicted survival probabilities for all training subjects are in object$survival
            #  out-of-bag prediction in object$survival.oob
        if "IPCW.times" in self.what:
            pass
            # self.times = predictRisk(fit, newdata=wdata, times=times)
        else:
            self.times = None
        # weights at subject specific event times
        if "IPCW.subject.times" in self.what:
            # pmat = fit.survival
            # jtimes = fit.time.interest
            # self.subject.times = pd.apply(1:len(subject.times), lambda(i: Ci = subject.times[i]))
            self.subject_times = pd.DataFrame()

            # pos = prodlim.sindex({'jump.times': jtimes, 'eval.times': Ci, 'comp': "smaller", 'strict': (lag == 1)) c(1, pmat[i,])[1 + pos]})
            pos = []
        else:
            self.subject_times = None
            out = {'times': self.times, 'subject.times': self.subject_times, 'method': self.method}

        out = self.output(out, self.keep, self.times, self.fit, self.call)
        return out

    def forest(self):
        # call = match.call()  ## needed for refit in crossvalidation loop
        #EHF = prodlim.EventHistory.frame(self.formula,self.data,specials=None,unspecialsDesign=False)
        EHF = {}
        # wdata = data.frame(cbind(unclass(EHF.event.history), EHF.design))
        wdata = {}
        ## wdata$status <- 1-wdata$status
        # wform = update(formula, "Surv(time,status)~.")
        wform = {}
        ## require(randomForestSRC)
        # if not NROW(na.omit(wdata)) > 0):
        #    raise('Stop')

        if not self.args:
            args = {'ntree':1000, 'importance':None}
            # fit = do.call(randomForestSRC.rfsrc, c(list(wform, data=wdata), args))
            fit = {}
            ## print(fit)
            fit.call = None
            # forest weights
            # FW = predict(fit, newdata=wdata, forest.wt = TRUE).forest.wt
            FW = {}
            #  weights at requested times
            #  predicted survival probabilities for all training subjects are in object$survival
            #  out-of-bag prediction in object$survival.oob
        if "IPCW.times" in self.what:
            # reverse Kaplan-Meier with forest weights
            # self.times = apply(data, 1, lambda (i: predict(prodlim.prodlim(Hist(time, status)~1, data = wdata, reverse = TRUE, caseweights = FW[i,]), times = times)
            self.times = []
        else:
            IPCW.times = None
            #  weights at subject specific event times

        if "IPCW.subject.times" in self.what:
            # self.subject_times = sapply(1:length(subject.times), lambda (i: prodlim:: predictSurvIndividual(prodlim::prodlim(Hist(time, status)~1, data = wdata, reverse = TRUE, caseweights = FW[i,]), lag = 1)[i])
            self.subject_times = []
        else:
            self.subject.times = None

        out = {'times': self.times,
               'subject_times': self.subject_times,
               'method': self.method}

        out = self.output(out, self.keep, self.times, self.fit, self.call)

        return out

    # reverse Kaplan-Meier, this is the one used in Metrics so focus on this.
    def marginal(self):
        # self.formula = update.formula(self.formula, "~1")
        self.formula = self.formula

        #reverse KM estimator (Schemper and Smith, 1996), that is the KM method with the
        # event indicator reversed so that the outcome of interest becomes being censored.
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2394262/
        # http://staff.pubhealth.ku.dk/~tag/Teaching/share/R-tutorials/SurvivalAnalysis.html
        # An approximation of this so-called reverse Kaplan-Meier can be obtained by reversing the event indicator
        # so that the outcome of interest becomes being censored.

        

        self.kmf.fit(self.data['T'].values, event_observed=self.data['failure_time'])
        self.kmf.survival_function_
        self.kmf.median_
        # fit = prodlim.prodlim(self.formula, data=data, reverse=True)
        fit = {}


        #  weights at requested times
        if "IPCW.times" in self.what:
            # self.times = predict(fit, newdata=data, times=times, level_chaos=1, mode="matrix", type="surv")
            self.times = []
        else:
            self.times = None

        # weights at subject specific event times
        if "IPCW.subject.times" in what:
            # self.subject_times = prodlim.predictSurvIndividual(fit, lag=self.lag)
            self.subject_times = []
        else:
            self.subject_times = None

        out = {'times': self.times,
               'subject_times': self.subject_times,
               'method': self.method}
        out = self.output(out, self.keep, self.times, self.fit, self.call)

        # class(out) < - "IPCW"
        return out
        ##   locsubject.times <- match(subject.times,fit$time,nomatch=NA)
        ##   if (any(is.na(locsubject.times))) stop("Can not locate all individual observation times" )
        ##   IPCW.subject.times <- c(1,fit$surv)[locsubject.times] ## at (subject.times_i-)
        ##   IPCW.times <- c(1,fit$surv)[prodlim::sindex(jump.times=fit$time,eval.times=times) +1] ## at all requested times

    # reverse Stone-Beran
    def nonpar(self):
        # call < - match.call()
        # fit = prodlim.prodlim(formula, data=data, reverse=TRUE, bandwidth="smooth")
        fit = {}
        #  weights at requested times
        if "IPCW.times" in self.what:
            # self.times = predict(fit, newdata=self.data, times=self.times, level.chaos = 1, mode = "matrix", type = "surv")
            self.times = []
        else:
            self.times = None
        # weights at subject specific event times
        if "IPCW.subject.times" in self.what:
            # self.subject_times = prodlim.predictSurvIndividual(fit, lag=self.lag)
            self.subject_times = []
        else:
            self.subject_times = None

        out = {'times': self.times,
               'subject_times': self.subject_times,
               'method': self.method}

        out = self.output(out, self.keep, self.times, self.fit, self.call)
        # class(out) < - "IPCW"
        return out

    # reverse Cox via Harrel's package
    def cox(self):
        # call < - match.call()
        # EHF = prodlim.EventHistory.frame(formula, data, specials=c("strat"), stripSpecials=c("strat"), specialsDesign=False,unspecialsDesign=False)
        EHF = {}
        if not EHF['strat']:
            # wdata = data.frame(cbind(unclass(EHF.event.history), EHF.design))
            wdata = {}
        else:
            wdata = {}
            # wdata = data.frame(cbind(unclass(EHF.event.history), EHF.design, EHF.strat))
        ## wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
        wdata.status = 1 - wdata.status
        # wform = update(formula, "Surv(time,status)~.")
        wform = {}
        # if not NROW(na.omit(wdata)) > 0):
        #    raise("Error")
        if not args:
            args = {'x': True, 'y': True, 'eps': 0.000001, 'surv': True}

        # fit = do.call(rms.cph, c(list(wform, data=wdata), args))
        fit = {}
        ## fit <- rms::cph(wform,data=wdata,surv=TRUE,x=TRUE,y=TRUE)

        #  weights at requested times
        if "IPCW.times" in self.what:
            # self.times = rms.survest(fit, newdata=wdata, times=times, se.fit = FALSE).surv
            self.times = []
        else:
            self.times = None

        # weights at subject specific event times
        if "IPCW.subject.times" in self.what:
            if self.lag == 1:
                # self.subject_times = rms.survest(fit,times=self.subject_times - min(diff(c(0, unique(self.subject_times)))) / 2, what='parallel')
                self.subject_times = []
            elif self.lag == 0:
                # self.subject_times = rms.survest(fit, times=self.subject_times, what='parallel')
                self.subject_times = []
            else:
                raise ("SubjectTimesLag must be 0 or 1")

        else:
            self.subject_times = None

        out = {'times': self.times, 'subject_times': self.subject_times, 'method': self.method}

        out = self.output(out, self.keep, self.times, self.fit, self.call)

        # class(out) < - "IPCW"
        return out

        # reverse Aalen method via the timereg package
        ## ipcw.aalen <- function(formula,data,method,args,times,subject.times,lag,what,keep){
        ## if (missing(lag)) lag=1
        ## if (missing(what)) what=c("IPCW.times","IPCW.subject.times")
        ## call <- match.call()
        ## EHF <- prodlim::EventHistory.frame(formula,
        ## data,
        ## specials=NULL,
        ## unspecialsDesign=FALSE)
        ## wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
        ## ## wdata <- as.data.frame(EHF)
        ## wdata$status <- 1-wdata$status
        ## wform <- update(formula,"Surv(time,status)~.")
        ## stopifnot(NROW(na.omit(wdata))>0)
        ## fit <- do.call(timereg::aalen,list(formula=formula,data=wdata,n.sim=0))
        ## fit$call <- NULL
        ## #  weigths at requested times
        ## if (match("IPCW.times",what,nomatch=FALSE)){
        ## IPCW.times <- predictRisk(fit,newdata=wdata,times=times)
        ## }  else {
        ## IPCW.times <- NULL
        ## }
        ## if (match("IPCW.subject.times",what,nomatch=FALSE)){
        ## if (lag==1)
        ## IPCW.subject.times <- diag(predictRisk(fit,newdata=data,times=pmax(0,subject.times-min(diff(unique(subject.times)))/2)))
        ## else if (lag==0)
        ## IPCW.subject.times <- diag(predictRisk(fit,newdata=data,times=subject.times))
        ## else stop("SubjectTimesLag must be 0 or 1")
        ## }
        ## else
        ## IPCW.subject.times <- NULL
        ## out <- list(times=times,IPCW.times=IPCW.times,IPCW.subject.times=IPCW.subject.times,fit=fit,call=call,method=method)
        ## class(out) <- "IPCW"
        ## out
        ## }
