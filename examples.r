#source("./metrics.r")

## Testing full results
calculate <- dget("metrics.r")
test = read.csv("test_data.csv", header = TRUE)
test['predictionscore'] <- order(test['D_PFS'])
pred = read.csv("predictions.csv", header = TRUE)
calculate(pred, test['D_PFS'], test['D_PFS_FLAG'], pred['study'])


## Testing prodLim package kaplan-meier reverse

