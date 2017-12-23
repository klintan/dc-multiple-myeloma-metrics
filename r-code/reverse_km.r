library(pec)
data(GBSG2,package="pec")
GBSG2$tgrade <- as.factor(as.character(GBSG2$tgrade))
GBSG2$Age <- cut(GBSG2$age,c(0,40,60,81),labels=c("21-40","41-60","61-80"))

library(prodlim)
library(survival)

quantile(prodlim(Hist(time,cens)~1,data=GBSG2,reverse=TRUE))
