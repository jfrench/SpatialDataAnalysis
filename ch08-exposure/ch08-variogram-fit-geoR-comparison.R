# fit a variogram model manually to a variogram

library(geoR)

# load smoky dataframe
load("./data/smoky.rda")

# create geodata object for geoR
geosmoky = as.geodata(cbind(smoky$longitude, smoky$latitude, smoky$ph))

# estimated semivariogram
vhat = variog(geosmoky, max.dist = 90, uvec = 20)

# extract relevant aspects of vhat
h <- vhat$u # distances
gammahat <- vhat$v # gammahat
nbins <- vhat$n # number of pairs in each bin

# fit exponential model using geoR package
fitexp = variofit(vhat, ini.cov.pars = c(.25, 30),
                  nugget = 0.05,
                  cov.model = "exponential",
                  weights = "cressie")
fitexp$cov.pars # estimated c = 0.283, a = 86.94
fitexp$nugget #estimated c0 = 0.085
fitexp$value/2 # the wrss

# create objective function
# theta = c(c, a, c0)
vfit <- function(theta) {
  # evaluate exponential model
  gammamod <- theta[3] + theta[1] * (1 - exp(-h/theta[2]))
  # compute wrss
  sum(nbins/gammamod^2 * (gammahat - gammamod)^2)/2
}

# fit exponential model
optim(par = c(0.25, 30, 0.05),
      fn = vfit,
      lower = c(0.001, 0.001, 0.001),
      upper = c(1, 1000, 1),
      method = "L-BFGS-B")

