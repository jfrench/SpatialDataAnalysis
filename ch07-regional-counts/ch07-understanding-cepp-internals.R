# understanding cepp.test internals
library(smerc)
# define information needed in function
data(nydf)
coords = as.matrix(nydf[,c("x", "y")])
cases = nydf$cases
pop = nydf$population
# expected number of cases in each region
ex = sum(cases)/sum(pop) * pop
# persons at risk
nstar = 10000
# don't use great circle distance
longlat = FALSE
# number of simulated data sets
nsim = 49
# simulation approach
simdist = "multinomial"

# taken from internals of cepp.test

# compute intercentroid distance
d = gedist(coords, longlat = longlat)
# determine the number of nearest neighbors to get nstar persons at risk
nn = casewin(d, pop, nstar)
# e.g., starting at region 1, you need to include regions 1, 2, and 15 to get
# nstar cases
nn[[1]]
# how much population from each region to add to get to nstar cases
# the population of a whole region is added before population from a new region
# is added
wts = cepp.weights(nn, pop, nstar)
# add proportions of cases for each candidate zone
cstar = sapply(seq_along(nn), function(i) {
  sum(cases[nn[[i]]] * wts[[i]])
}, USE.NAMES = FALSE)
# do the same thing for simulated data sets
# returns the largest cstar for each simulated data set
csim = cepp.sim(nsim = nsim, nn = nn, ty = sum(cases), ex = ex,
                wts = wts, simdist = simdist)
# compute p-value for each candidate zone
pvalue = mc.pvalue(cstar, csim)

# most likely cluster
which.max(cstar)

# rank clusters by largest test statistic
order(cstar, decreasing = TRUE)

# view statistics of each window by ranking
cstar[order(cstar, decreasing = TRUE)]
