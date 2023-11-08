# fit a variogram model manually to a variogram

library(gstat)

# turn smoky into sf dataframe
load("./data/smoky.rda")
smoky <- sf::st_as_sf(smoky,
                      coords = c("longitude", "latitude"))
# create gstat object for further analysis
# formula ph ~ 1 assumes a constant mean over the spatial domain
gsmoky = gstat(id = "ph", formula = ph ~ 1, data = smoky)

# compute empirical semivariogram
vhat = variogram(gsmoky, cutoff = 90, width = 4.5)
plot(vhat)


# extract relevant aspects of vhat
h <- vhat$dist # distances
gammahat <- vhat$gamma # gammahat
nbins <- vhat$np # number of pairs in each bin

### Table 8.1 fit variogram model using WRSS
# estimated exponential model with starting values c = .25, a = 30, c0 = .05
fitexp = fit.variogram(vhat,
                       vgm(.25, "Exp", 30, .05),
                       fit.method = 2)
fitexp #c = 0.211, a = 33.18, c0 = 0.042
# plot variogram with estimated exponential model
plot(vhat, fitexp, main = "WRSS exponential fit")
# wrss of fit
attr(fitexp, "SSErr")

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

