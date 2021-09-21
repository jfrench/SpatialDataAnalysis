# install.packages(spdep", "sf")
library(spdep)
library(sf)

# read shapefile for new york counties
ny8 <- sf::st_read("./data/NY_data/ny8_utm18.shp")
# read neighbor information
ny_nb <- spdep::read.gal("./data/NY_data/NY_nb.gal", override.id = TRUE)

# plot region boundaries from ny8
plot(st_geometry(ny8), border="grey60")
# plot neighbors
plot(ny_nb, coords = st_centroid(st_geometry(ny8)),
     add=TRUE, col="blue", pch = 19, cex = 0.6)

### geary's c
# assume adjacency weights (w_ij = 1 if regions i and j share a boundary)
# proximity matrix, binary style.  W is row standardized.
w = nb2mat(ny_nb, style = "B")
# see ?nb2listw for more options
# proximaty matrix in list format
lw = nb2listw(ny_nb, style = "B")

geary.test(ny8$Cases, listw = lw, randomisation = FALSE)
# base test w/ randomization p-value
geary.mc(ny8$Cases, listw = lw, nsim = 499)

# base test w/ Monto Carlo p-value, simulating data under constant risk hypothesis
# some preliminaries
N = length(ny8$Cases) # number of regions
y = ny8$Cases # number of cases
n = ny8$POP8 #population sizes
r <- sum(y)/sum(n) # estimated risk
rni <- r * n # expected per region

# observed geary's statistic
nsim = 499
t0 = geary(y, listw = lw, n = N, n1 = N - 1, S0 = Szero(lw))$C
# simulate data under CRH
tsim = numeric(nsim)
# calculate geary's c for poisson data simulated under crh
for (i in 1:nsim) {
  tsim[i] = geary(rpois(N, rni), listw = lw, n = N, n1 = N - 1, S0 = Szero(lw))$C
}

# p-value for geary's c constant risk monte carlo test
(sum(tsim <= t0) + 1)/(nsim + 1)

## Do the same tests with incidence rates
rates = y/n
# incidence test w/ normality approximation for p-value
geary.test(rates, listw = lw, randomisation = FALSE)
# incicdence test w/ randomization p-value
geary.mc(rates, listw = lw, nsim = 499)

# incidence test w/ Monto Carlo p-value, simulating data under constant risk hypothesis
# some preliminaries
# observed geary's statistic
t0b = geary(rates, listw = lw, n = N, n1 = N - 1, S0 = Szero(lw))$C
# calculate geary's c for poisson data simulated under crh, after rate
tsimb = numeric(nsim)
# correction
for (i in 1:nsim) {
  tsimb[i] = geary(rpois(N, rni)/n, listw = lw, n = N, n1 = N - 1, S0 = Szero(lw))$C
}

# p-value for geary's c constant risk monte carlo test for incidence rate
(sum(tsimb <= t0b) + 1)/(nsim + 1)

#### summary of results
### Counts
## normality assumption: conclude spatial autocorrelation
## randomization assumption: conclude spatial autocorrelation
## Monte Carlo CRH assumption:  no spatial autocorrelation
# Conclusion: similarity in counts caused by heterogeneities in population
# size, not similarities in spatial deviation from the mean

### Incidence proportion
## normality assumption: conclude spatial autocorrelation
## randomization assumption: no spatial autocorrelation
## Monte Carlo CRH assumption:  suggestive, not conclusive autocorrelation

## results are sensitive to assumptions.
