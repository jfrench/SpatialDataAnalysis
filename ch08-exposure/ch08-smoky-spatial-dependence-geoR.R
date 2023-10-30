# install.packages(c(""geoR", "autoimage"))
library(geoR) # estimating anisotropic covariance models

### some notes
# the variog and likfit functions (from the geoR package)
# required geodat objects, which are created using the
# as.geodata function.

# the variogram functions for the gstat package parameterize
# the direction in degrees.  the variogram functions for the
# geoR package function the variogram functions in radians,
# where radians = degrees/180*pi

# load smoky dataframe
load("./data/smoky.rda")

# create geodata object for geoR
geosmoky = as.geodata(cbind(smoky$longitude, smoky$latitude, smoky$ph))

# bubble plot with scale
autoimage::autopoints(smoky$longitude, smoky$latitude, smoky$ph,
                      xlab = "longitude", ylab = "latitude")

### Fig 8.5
# construct variogram up to distance of 90 with 20 bins (width = 4.5)
variog1b = variog(geosmoky, max.dist = 90, uvec = 20)
plot(variog1b)

### Fig 8.6
# construct variogram up to distance of 120 with 12 bins (widths = 10)
variog2b = variog(geosmoky, max.dist = 120, uvec = 12)
plot(variog2b)

# number of pairs used to estimate semivariance in each interval
variog2b$n # geoR

### Table 8.1 fit variogram model using WRSS
# estimate exponential model with starting values c = .25, a = 30, c0 = .05
# ini.cov.pars = c(c, a). nugget = c0
fitexpb = variofit(variog2b, ini.cov.pars = c(.25, 30),
                   nugget = 0.05, cov.model = "exponential",
                   weights = "cressie")
fitexpb # estimated c0 = 0.0277, c = 0.2070, a = 25.1576
fitexpb$value # wrss of fit

# plot empirical semivariogram
plot(variog2b)
# add fitted model to empirical semivariogram
lines(fitexpb)

### Table 8.2 # fit variogram model using WRSS
# estimated spherical model with starting values c = .25, a = 30, c0 = .05
fitsphb = variofit(variog2b, ini.cov.pars = c(.25, 30), nugget = 0.05,
                   cov.model = "spherical", weights = "cressie")
fitsphb  # estimated c0 = 0.0558, c = 0.1683, a = 62.6156
fitsphb$value # wrss of fit
# the spherical model fits worse since its WRSS = 30.44
# vs 30.06 for the exponential model
# plot empirical semivariogram
plot(variog2b)
# add fitted model to empirical semivariogram
lines(fitsphb)

# estimated matern model with starting values c = .25, a = 30, c0 = .05
fitmatb = variofit(variog2b, ini.cov.pars = c(.25, 30),
                   nugget = 0.05,
                   kappa = 0.75,
                   cov.model = "matern",
                   weights = "cressie",
                   fix.kappa = FALSE)
fitmatb  # estimated c0 = 0.0529, c = 0.1761, a = 16.36, smoothness =
fitmatb$value # wrss of fit 29.875
# the matern model fits the best since it has the smallest RSS
plot(variog2b)
# add fitted model to empirical semivariogram
lines(fitmatb)

# fit covariance model using Restricted Maximum Likelihood

# REML estimation covariance parameters of geosmoky for an exponential model
# with starting values c = .25, a = 30, c0 = 0.05
lfit_exp = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05,
                  cov.model = "exponential", lik.method = "REML")
lfit_exp # c = .1929, a = 10.97, c0 = 0, muhat = 7.1929

# estimate an isotropic matern covariance model with and estimating the
# smoothness parameter (kappa). Note fix.kappa = FALSE.'
# with starting values c = .25, a = 30, c0 = 0.05, smoothness = 0.5
lfit_mat = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05,
                  cov.model = "matern", lik.method = "REML", kappa = 1,
                  fix.kappa = FALSE)
lfit_mat # c = .1841, a = 7.1482, c0 = 0, smoothness = 0.6872, muhat = 7.1663

### Fig 8.11
# plot directional variogram
# direction are 70, 115, 160, 205 (205 - 180 = 25) degrees
# same (general) thing using geoR.  Angles must be between 0 and 180.  I get 25 from 205 - 180
# geoR specifies direction in radians, which is degrees/180*pi
variog3b = variog4(geosmoky, direction = c(25, 70, 115, 160)/180*pi, max.dist = 120)
plot(variog3b)

# geoR can't estimate anisotropic semivariogram models, so we must use
# likelihood-based methods

# fit geometric anisotropy covariance model using REML
# psiA controls the anisotropy angle. This is the direction of amin.
# psiR controls the ratio of the large a (amaj) to the smallest a (amin)
# note that psiA must be converted from degrees to radians using the formula radians = degrees/180*pi
# by default psiA and psiR are considered FIXED (not estimated)
# in model 1, both are fixed, in model 2, the psiA is fixed, and in model 3, neither are fixed

# estimate an anisotropic exponential covariance model with fixed angle and ratio
# 70 is the angle chosen in the book example.t
# the /180*pi converts the angle to radians
# ini.cov.pars = c(c, amaj), nugget = c0
# initial angle of 70 degrees, initial amaj/amin = 3
lfit1 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05,
               cov.model = "exponential", lik.method = "REML",
               psiA = 70/180*pi, psiR = 3)
lfit1 # muhat = 7.2349, c0 = 0.0044, c = 0.2198, amin = 9.4746, psiA = 70 degrees, psiR = 3
# estimate an anisotropic exponential covariance model with fixed angle but also estimates
# the ratio of amaj/amin (amax/amin). Use fix.psiR = FALSE.
lfit2 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05,
               cov.model = "exponential", lik.method = "REML",
               psiA = 70/180*pi, psiR = 3, fix.psiR = FALSE)
lfit2 # muhat = 7.19, c0 = 0.0016, c = 0.1880, amin = 10.82, psiA = 70 degrees, psiR = 1.91

# estimate an anisotropic exponential covariance model with both angle and ratio
# being estimated. Note fix.psiR (ratio) = FALSE, and fix.psiA = FALSE.
lfit3 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05,
               cov.model = "exponential", lik.method = "REML",
               psiA = 70/180*pi, psiR = 3, fix.psiR = FALSE, fix.psiA = FALSE)
lfit3 # muhat = 7.21, c0 = 0.0039, c = 0.2013, amin = 12.83, psiA = 0.92 radians, psiR = 3.07

# estimate an anisotropic matern covariance model with smoothness, angle, and ratio
# being estimated. Note fix.kappa = FALSE (smoothness).
# Note fix.psiR (ratio) = FALSE, and fix.psiA = FALSE.
lfit4 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05, kappa = 1,
               cov.model = "matern", lik.method = "REML",
               psiA = 70/180*pi, psiR = 3, fix.psiR = FALSE, fix.psiA = FALSE,
               fix.kappa = FALSE)
lfit4 # muhat = 7.61, c0 = 0, c = 2.54, amin = 903.01, kappa = 0.35, psiA = 0.85 radians, psiR = 2.93

# compare AIC/BIC values.  Model 3 is best in terms of AIC since it has the
# smallest AIC, but we'll use model 2 because that's what the book uses
summary(lfit1) # AIC = 60.92
summary(lfit2) # AIC = 61.81, c = .188, amin = 10.82, c0= 0.0016, amaj/amin = 1.91, degrees = 70
summary(lfit3) # AIC = 60.73
summary(lfit4) # AIC = 62.29