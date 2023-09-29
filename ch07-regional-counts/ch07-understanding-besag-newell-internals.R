# understanding bn.test internals
library(smerc)

# perform besag-newell test for c* = 23 cases
cstar <- 23
# significance level
alpha <- 0.05
# approach for computing pvalue
modified <- FALSE

# define information needed in function
data(nydf)
coords <- as.matrix(nydf[,c("x", "y")])

cases <- nydf$cases
pop <- nydf$population
# expected number of cases in each region
ex <- sum(cases)/sum(pop) * pop

# get distances
d <- gedist(coords)
# compute case windows
cwins <- casewin(d, cases, cstar)
# determine number of retions in each case window
l <- sapply(cwins, length)
# determine number of cases in each case window
case_cwins <- zones.sum(cwins, cases)
# determined expected number of cases in each case window
ex_cwins <- zones.sum(cwins, ex)

# compute pvalue depending on whether the modified test is
# being perform
if (!modified) {
  pvalue <- stats::ppois(cstar - 1, lambda = ex_cwins,
                         lower.tail = FALSE)
} else {
  pvalue <- stats::ppois(case_cwins - 1, lambda = ex_cwins,
                         lower.tail = FALSE)
}

# determine most significant non-overlapping clusters
pruned <- sig_noc(tobs = l, zones = cwins, pvalue = pvalue,
                  alpha = alpha, order_by = "pvalue")

pruned
