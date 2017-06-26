#######################################################################################
# options, libraries

options(stringsAsFactors=FALSE)

library(plotly)
library(broom)

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
# check that estimated OR checks out 
pCases <- .75
pControls <- .25
pCases/(1-pCases)/(pControls/(1-pControls))
summary(unlist(lapply(dats, function(D) {
   exp(coef(glm(y~x,family=binomial, data=D)))[2]})))

# clean up 
rm(list=setdiff(ls(), "binSim"))

#######################################################################################
# run simulation
# plan: fix n, pCases, and cacoRatio in a few reasonable panels; then,
# plot estimate CIs over range of pControls

iter <- 250
pControls <- seq(.05,.95,by=.02)
sim <- list()

### pCases=.95
n <- 300
cacoRatio <- 2
pCases <- .95
set.seed(8675309)
for (k in 1:length(pControls)) {
   sim[[k]] <- data.frame(est=rep(0, iter), lower=0, upper=0, 
                          OR=0, ORlower=0, ORupper=0, p=0)
   for (i in 1:iter) {
      dat <- binSim(n, cacoRatio, pCases, pControls=pControls[k])
      fit <- glm(y ~ x, family=binomial, data=dat)
      tidyobj <- tidy(fit)[2,]
      tidyCI <- confint_tidy(fit)[2,]
      sim[[k]][i,] <- c(tidyobj$est, tidyCI,
         exp(tidyobj$est), exp(tidyCI), tidyobj$p.value)
   }
}

simdat <- data.frame(est=unlist(lapply(sim, function(x) median(x$est))),
                      lower=unlist(lapply(sim, function(x) median(x$lower))),
                      upper=unlist(lapply(sim, function(x) median(x$upper))),
                      OR=unlist(lapply(sim, function(x) median(x$OR))),
                      ORlower=unlist(lapply(sim, function(x) median(x$ORlower))),
                      ORupper=unlist(lapply(sim, function(x) median(x$ORupper))),
                      ps=unlist(lapply(sim, function(x) mean(x$p))),
                      pControls=pControls, pCases=pCases, 
                      ORtrue=pCases/(1-pCases)/(pControls/(1-pControls)))



### pCases=.75
n <- 300
cacoRatio <- 2
pCases <- .75
set.seed(8675309)
for (k in 1:length(pControls)) {
  sim[[k]] <- data.frame(OR=rep(0, iter), lower=0, upper=0, p=0)
  for (i in 1:iter) {
    dat <- binSim(n, cacoRatio, pCases, pControls=pControls)
    fit <- glm(y ~ x, family=binomial, data=dat)
    tidyobj <- tidy(fit)[2,]
    tidyCI <- confint_tidy(fit)[2,]
    sim[[k]][i,] <- c(exp(tidyobj$est), exp(tidyCI), tidyobj$p.value)
  }
}

# bind this simulated data to the previous
binddat <- data.frame(OR=unlist(lapply(sim, function(x) mean(x$OR))),
                     lower=unlist(lapply(sim, function(x) mean(x$lower))),
                     upper=unlist(lapply(sim, function(x) mean(x$upper))),
                     ps=unlist(lapply(sim, function(x) mean(x$p))),
                     pControls=pControls, pCases=pCases)
simdat <- rbind(simdat, binddat)

### pCases=.5
n <- 300
cacoRatio <- 2
pCases <- .5
set.seed(8675309)
for (k in 1:length(pControls)) {
  sim[[k]] <- data.frame(OR=rep(0, iter), lower=0, upper=0, p=0)
  for (i in 1:iter) {
    dat <- binSim(n, cacoRatio, pCases, pControls=pControls[k])
    fit <- glm(y ~ x, family=binomial, data=dat)
    tidyobj <- tidy(fit)[2,]
    tidyCI <- confint_tidy(fit)[2,]
    sim[[k]][i,] <- c(exp(tidyobj$est), exp(tidyCI), tidyobj$p.value)
  }
}

# bind this simulated data to the previous
binddat <- data.frame(OR=unlist(lapply(sim, function(x) mean(x$OR))),
                      lower=unlist(lapply(sim, function(x) mean(x$lower))),
                      upper=unlist(lapply(sim, function(x) mean(x$upper))),
                      ps=unlist(lapply(sim, function(x) mean(x$p))),
                      pControls=pControls, pCases=pCases)
simdat <- rbind(simdat, binddat)

#ggplot2 can't handle infinite or very large numbers
simdat$upper[which(is.infinite(simdat$upper))] <- 1e10
simdat$upper[which(is.na(simdat$upper))] <- 1e10
simdat$upper[which(simdat$upper>1e10)] <- 1e10 

ORp <- ggplot(simdat, aes(y=OR, x=pControls, colour=as.factor(pCases))) + 
   geom_smooth(se=FALSE) + 
   geom_ribbon(aes(ymin=lower, ymax=upper, fill="band"), alpha = 0.3) +
   xlab("P(Exposed | Control)") + ylab("Estimated odds ratio") +
   scale_y_continuous(limits=c(0, 5))
   #coord_cartesian(ylim=c(0, 10)) 
   #facet_wrap(~pCases) #invokes side-by-side panels

ggplotly(ORp)

simdat$ORupper[which(is.na(simdat$ORupper))] <- 1e10
plot(simdat$pControls, simdat$OR, ylim=range(c(simdat$ORlower, simdat$ORupper)))
lines(simdat$pControls, simdat$ORlower, lty=2)
lines(simdat$pControls, simdat$ORupper, lty=2)

plot(diffs, simdat$OR)
