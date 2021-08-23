# install.packages(spdep", "rgdal",)
library(spdep)
library(rgdal)

# read in data
setwd("~/OneDrive - The University of Colorado Denver/Teaching/Math4027/Data/NY_data")
NY8 <- readOGR(".", "NY8_utm18")
NY_nb <- read.gal("NY_nb.gal", region.id = row.names(NY8))

setwd("~/OneDrive - The University of Colorado Denver/Teaching/Math4027/Data/")
nydf = read.table("NYtract.dat")
names(nydf) = c("x", "y", "Population", "Observed")

### geary's c
# assume adjacency weights (w_ij = 1 if regions i and j share a boundary)
# proximity matrix, binary style.  W is row standardized.
W = nb2mat(NY_nb, style = "B")
# see ?nb2listw for more options
# proximaty matrix in list format
lw = nb2listw(NY_nb, style = "B")

# base test w/ normality approximation for p-value
geary.test(NY8$Cases, listw = lw, randomisation = FALSE)
# base test w/ randomization p-value
geary.mc(NY8$Cases, listw = lw, nsim = 499)

# base test w/ Monto Carlo p-value, simulating data under constant risk hypothesis
# some preliminaries
N = length(NY8$Cases) # number of regions
y = NY8$Cases # number of cases
n = NY8$POP8 #population sizes
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
