# install.packages(c("spdep", "rgdal", "RColorBrewer", "colorRamps"))
library(spdep)
library(sf)
library(spatialreg)
library(ggplot2)
library(patchwork)

# read in data
# setwd("~/.data/NY_data")
NY8 <- sf::st_read("./data/NY_data/NY8_utm18.shp")
TCE <- sf::st_read("./data/NY_data/TCE.shp")
NY_nb <- spdep::read.gal("./data/NY_data/NY_nb.gal", region.id = row.names(NY8),
                         override.id = TRUE)
cities <- sf::st_read("./data/NY_data/NY8cities.shp")

# plot of regions with TCE sites (and name)
par(mfrow = c(1, 2))
plot(NY8$geometry, border="grey60", axes = TRUE, main = "cities")
text(st_coordinates(cities), labels=as.character(cities$names), font=2, cex=0.6)
plot(NY8$geometry, border="grey60", axes=TRUE, main = "factories")
plot(TCE$geometry, pch=19, cex=0.7, add = TRUE)
text(st_coordinates(TCE), labels=as.character(TCE$name), cex=0.7,
     font=1, pos=c(4,1,4,1,4,4,4,2,3,4,2), offset=0.3)
par(mfrow = c(1, 1))

# fit ols lm based on 3 covariates
nylm <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = NY8)
summary(nylm) # exposure to TCE doesn't appear important, though the
# other census variables are.

# add residuals from lm fit to NY8
NY8$lmresid <- residuals(nylm)

# obtain binary spatial proximity matrix
# in proper format for NY8 data
NYlistw <- nb2listw(NY_nb, style = "B")

# some evidence of correlated errors, i.e., residual dependence
lm.morantest(nylm, NYlistw)

# fit one-parameter SAR model with same covariates
nysar <- spatialreg::spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                  data = NY8, listw = NYlistw)
summary(nysar)
# conclude that (lambda, i.e., rho in our book) is different
# from zero according to a LRT
# exposure to TCE still not significant, but close enough
# that it shouldn't be discarded
# other two variables associated with higher incidence rates

# plot trend and stochastic components
# non-spatial component of fitted trend
NY8$sar_trend <- nysar$fit$signal_trend
# spatial component of fitted trend
NY8$sar_stochastic <- nysar$fit$signal_stochastic

# fill colors the regions according to the specified variable
# scale_fill_gradient2 uses diverging color schole with a white midpoint at 0
# limits is not usually needed. It is included so that the color scale
# is the same across the various plots below.
a <- ggplot(NY8) + geom_sf(aes(fill = sar_trend)) +
  scale_fill_gradient2(limits = c(-1, 1.3))
b <- ggplot(NY8) + geom_sf(aes(fill = sar_stochastic)) +
  scale_fill_gradient2(limits = c(-0.16, 0.39))
# use patchwork to plot a and b side-by-side
a + b

# changes weights from ls regression.  Weights are inversely proportional
# to population sizes in each region.
nylmw <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = NY8,
            weights = POP8)
summary(nylmw)
# TCE exposure significant with expected sign!

# add residual from wls fit to NY8
NY8$lmwresid <- residuals(nylmw)

# compare residuals from previous 2 models
a2 <- ggplot(NY8) + geom_sf(aes(fill = lmresid)) +
  scale_fill_gradient2()
b2 <- ggplot(NY8) + geom_sf(aes(fill = lmwresid)) +
  scale_fill_gradient2()
# use patchwork to plot a and b side-by-side
a2 + b2

# check moran's after accounting for heterogeneous population sizes
lm.morantest(nylmw, NYlistw)
# model mispecification (rho != 0) in first test caused by
# heteroskedasticity more than to spatial autocorrelation?

# SAR model with different weights, that account for heterogeneous
# population sizes.  Fit sar model with different weights for Vv.  Specifically,
# Vv propto 1/ni for each region i
nysarw <- spatialreg::spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                   data = NY8, listw = NYlistw, weights = POP8)
summary(nysarw)
# no traces of spatial autocorrelation after adjusting for the
# heterogeneous population sizes

# separate trend from stochastic component
NY8$sarw_trend <- nysarw$fit$signal_trend
NY8$sarw_stochastic <- nysarw$fit$signal_stochastic

a3 <- ggplot(NY8) + geom_sf(aes(fill = sarw_trend)) +
  scale_fill_gradient2(limits = c(-1, 1.3))
b3 <- ggplot(NY8) + geom_sf(aes(fill = sarw_stochastic)) +
  scale_fill_gradient2(limits = c(-0.16, 0.39))
# use patchwork to plot a and b side-by-side
a3 + b3
# residual autocorrelation vanishes

# fit one-parameter car model
nycar <- spatialreg::spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                  data = NY8, family = "CAR", listw = NYlistw)
summary(nycar)
# the rho parameter is significant.  Suggests spatial
# autocorrelation

# fit car model with different weights for Vc.  Specifically,
# Vc propto 1/ni for each region i
nycarw <- spatialreg::spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                   data = NY8, family = "CAR", listw = NYlistw, weights = POP8)
summary(nycarw)
# spatial autocorrelation vanishes again

# extract fitted values from sarw and carw models
NY8$fitsarw = nysarw$fit$fitted.values
NY8$fitcarw = nycarw$fit$fitted.values

sar_car_limits <- range(c(NY8$fitsarw, NY8$fitcarw))

a4 <- ggplot(NY8) + geom_sf(aes(fill = fitsarw)) +
  scale_fill_viridis_c(limits = sar_car_limits)
b4 <- ggplot(NY8) + geom_sf(aes(fill = fitcarw)) +
  scale_fill_viridis_c(limits = sar_car_limits)
# use patchwork to plot a and b side-by-side
a4 + b4 # plot fitted values from sarw and carw (using two color palettes)

# compare fitted values, manually
NY8$fitsarw - NY8$fitcarw

# overall conclusion
# we probably didn't need to account for spatial autocorrelation, assuming
# we accounted for varying population sizes.  Spatial patterns driven by
# population sizes, not true residual spatial autocorrelation.
