#source("./metrics.r")

## Testing full results
#calculate <- dget("metrics.r")
test = read.csv("test_data.csv", header = TRUE)
test['predictionscore'] <- order(test['D_PFS'])
#pred = read.csv("predictions.csv", header = TRUE)
#calculate(pred, test['D_PFS'], test['D_PFS_FLAG'], pred['study'])


## Testing prodLim package kaplan-meier reverse

# T              : vector of observed failure times
# delta          : vector of indicator of status

library(prodlim)
library(rms)

times = 30.5 * c(14, 16, 18, 20, 22)

T = test['D_PFS']
colnames(T) <- c("T")

delta = test['D_PFS_FLAG']
colnames(delta) <- c("delta")

#T <- T[order(T$T),]


print(T)
order_T = order(T)
print(order_T)
T <- T[order_T,]
delta <- delta[order_T,]


# reverse the censored events (making them non-censored)
df = data.frame(failure_time=T,status=as.numeric(delta!=0))

# using the marginal Kaplan-Meier for the censoring times
weights <- pec::ipcw(Surv(failure_time,status)~1,data=df,method="marginal",times=times,subjectTimes=T,subjectTimesLag=1)

print(weights)