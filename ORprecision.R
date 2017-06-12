#######################################################################################
# options, libraries

options(stringsAsFactors=FALSE)

#######################################################################################

# function to simulate binomial data from which to calculate odds ratios
# n: overall sample size
# cacoRatio: average case-control ratio. Default assumes 2 controls per case
# pCases: expeted exposure prevalence in cases
# pControls: expected exposure prevalence in controls

binSim <- function(n, cacoRatio=2, pCases, pControls) {
   y <- x <- rep(0, n)
   nCases <- round(n/(cacoRatio+1))
   y[1:nCases] <- 1
   y[(nCases+1):n] <- 0
   x[y==1] <- rbinom(sum(y), 1, pCases)
   x[y==0] <- rbinom(n-sum(y), 1, pControls)
   return(data.frame(y=y, x=x))
}

# test
iter <- 1000
dats <- vector("list", length=iter)
for (i in 1:iter) {
   dats[[i]] <- binSim(n=250, cacoRatio=3, pCases=.75, pControls=.25)
}
# check ratio of controls:cases
summary(unlist(lapply(dats, function(D) {sum(D$y==0)/sum(D$y==1)})))
# check proportion of exposed in cases
summary(unlist(lapply(dats, function(D) {sum(D$x[D$y==1])/sum(D$y==1)})))
# check proportion of exposed in controls
summary(unlist(lapply(dats, function(D) {sum(D$x[D$y==0])/sum(D$y==0)})))
