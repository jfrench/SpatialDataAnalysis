# install.packages("spdep", "sf")
library(spdep)
library(sf)

# read shapefile for new york counties
ny8 <- sf::st_read("./data/NY_data/NY8_utm18.shp")
# read neighbor information
ny_nb <- spdep::read.gal("./data/NY_data/NY_nb.gal", override.id = TRUE)

# plot region boundaries from ny8
plot(st_geometry(ny8), border="grey60")
# plot neighbors
plot(ny_nb, coords = st_centroid(st_geometry(ny8)),
     add=TRUE, col="blue", pch = 19, cex = 0.6)

### moran's i
# assume adjacency weights (w_ij = 1 if regions i and j share a boundary)
# proximity matrix, binary style.  W is row standardized.
w = nb2mat(ny_nb, style = "B")
# see ?nb2listw for more options
# proximaty matrix in list format
lw = nb2listw(ny_nb, style = "B")

# base test w/ normality approximation for p-value
moran.test(ny8$Cases, listw = lw, randomisation = FALSE)
# base test w/ randomization p-value
(ir = moran.mc(ny8$Cases, listw = lw, nsim = 499))
# base test w/ Monto Carlo p-value, simulating data under constant risk hypothesis
# some preliminaries
N = length(ny8$Cases) # number of regions
y = ny8$Cases # number of cases
n = ny8$POP8 #population sizes
r <- sum(y)/sum(n) # estimated risk
rni <- r * n # expected per region

# observed moran's statistic
nsim = 499
t0 = moran(y, listw = lw, n = N, S0 = Szero(lw))$I
# simulate data under CRH
tsim = numeric(nsim)
# calculate moran's i for poisson data simulated under crh
for (i in 1:nsim) {
  tsim[i] = moran(rpois(N, rni), listw = lw, n = N, S0 = Szero(lw))$I
}

# p-value for moran's i constant risk monte carlo test
(sum(tsim >= t0) + 1)/(nsim + 1)

# compare histogram's for randomization and constant risk monte carlo test
hist(ir$res, xlab = "Moran's I", main = "Cases: Randomization")
abline(v = ir$statistic)

hist(tsim, xlab = "Moran's I", main = "Cases: Constant Risk")
abline(v = t0)

## Do the same tests with incidence rates
rates = y/n
# incidence test w/ normality approximation for p-value
moran.test(rates, listw = lw, randomisation = FALSE)
# incicdence test w/ randomization p-value
(irb = moran.mc(rates, listw = lw, nsim = 499))
# incidence test w/ Monto Carlo p-value, simulating data under constant risk hypothesis
# some preliminaries

# observed moran's statistic
t0b = moran(rates, listw = lw, n = N, S0 = Szero(lw))$I
# calculate moran's i for poisson data simulated under crh, after rate
tsimb = numeric(nsim)
# correction
for (i in 1:nsim) {
  tsimb[i] = moran(rpois(N, rni)/n, listw = lw, n = N, S0 = Szero(lw))$I
}

# p-value for moran's i constant risk monte carlo test for incidence rate
(sum(tsimb >= t0b) + 1)/(nsim + 1)

# compare histogram's for randomization and constant risk monte carlo test
# for incidence rates
hist(irb$res, xlab = "Moran's I", main = "Incidence Proportions: Randomization")
abline(v = irb$statistic)

hist(tsimb, xlab = "Moran's I", main = "Incidence Proportions: Constant Risk")
abline(v = t0b)

### Use CR Moran's I for inference
# make a function out of this process
i_cr = function(y, rni, w) {
  y_std = matrix((y - rni)/sqrt(rni))
  return(sum(w * y_std %*% t(y_std))/sum(w))
}

tsimc = numeric(nsim)
t0c = i_cr(y, rni, w) # observed statistic
# statistics for data simualted under CRH
for (i in 1:nsim) tsimc[i] = i_cr(rpois(N, rni), rni = rni, w = w)
# p-value
(sum(tsimc >= t0c) + 1)/(nsim + 1)

#### summary of results
### Counts
## normality assumption: conclude spatial autocorrelation
## randomization assumption: conclude spatial autocorrelation
## Monte Carlo CRH assumption:  no spatial autocorrelation
# Conclusion: similarity in counts caused by heterogeneities in population
# size, not similarities in spatial deviation from the mean

### Incidence proportion
## normality, randomization, Monte Carlo CRH assumption: no spatial autocorrelation
## this doesn't account for the fact that our rate estimates aren't as good
## for regions with smaller population sizes

### I_cr for counts
## Conclude clustering of standardized deviations of the
## regional counts
## from their expected values under the CRH

### Concluding thoughts
# 1. normality and randomization assumptions are inappropriate and
# can lead to wrong inference.
# 2.  Tests using counts under CRH and rates do not suggest
# clustering.
# 3.  Strong evidence for spatial correlation of
# standardized deviation from
# regional expected counts under CRH.
# 4.  The apparent discrepancies between 2 and 3 are caused the fact that
# we do not adjust for heterogeneous means/variances across space (the standard Moran's
# i using a common mean and standard deviation).  Accounting for this provides additional
# precision in assessing the spatial similarity of statistically unusual counts.
