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

#######################################################################################
# run simulation
# plan: fix n, pCases, and cacoRatio in a few reasonable panels; then,
# plot OR CI over range of pControls

iter <- 250
pControls <- seq(.05,.95,by=.02)
sim1 <- list()

n <- 300
cacoRatio <- 2
pCases <- .95
set.seed(8675309)
for (k in 1:length(pControls)) {
   sim1[[k]] <- data.frame(OR=rep(0, iter), lower=0, upper=0, p=0)
   for (i in 1:iter) {
      dat <- binSim(n, cacoRatio, pCases, pControls=pControls)
      fit <- glm(y ~ x, family=binomial, data=dat)
      tidyobj <- tidy(fit)[2,]
      tidyCI <- confint_tidy(fit)[2,]
      sim1[[k]][i,] <- c(exp(tidyobj$est), exp(tidyCI), tidyobj$p.value)
   }
}

sim1dat <- data.frame(OR=unlist(lapply(sim1, function(x) mean(x$OR))),
                      lower=unlist(lapply(sim1, function(x) mean(x$lower))),
                      upper=unlist(lapply(sim1, function(x) mean(x$upper))),
                      ps=unlist(lapply(sim1, function(x) mean(x$p))),
                      diffs=diffs)
sim1dat$upper[which(is.infinite(sim1dat$upper))] <- 1e10
sim1dat$upper[which(is.na(sim1dat$upper))] <- 1e10
sim1dat$upper[which(sim1dat$upper>1e10)] <- 1e10 #ggplot2 can't handle very large numbers

p <- ggplot(sim1dat) + geom_smooth(aes(y=OR, x=diffs), se=FALSE) +
   geom_ribbon(aes(ymin=lower, ymax=upper, x=diffs, fill = "band"), alpha = 0.3) +
   xlab("P(Exposed | Case) - P(Exposed | Control)") + ylab("Estimated odds ratio")
   #facet_wrap(~nGrp) #invokes side-by-side panels, if needed

ggplotly(p)


plot(diffs, sim1dat$OR)
