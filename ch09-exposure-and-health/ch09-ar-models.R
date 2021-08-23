# install.packages(c("spdep", "rgdal", "RColorBrewer", "colorRamps"))
library(spdep)
library(rgdal)
library(RColorBrewer)

# read in data
setwd("~/OneDrive - The University of Colorado Denver/Teaching/Math6384/Data/NY_data")
NY8 <- readOGR(".", "NY8_utm18")
TCE <- readOGR(".", "TCE")
NY_nb <- read.gal("NY_nb.gal", region.id=row.names(NY8))
cities <- readOGR(".", "NY8cities")

# plot of regions with TEC sites (and name)
plot(NY8, border="grey60", axes=TRUE)
text(coordinates(cities), labels=as.character(cities$names), font=2, cex=0.9)
text(bbox(NY8)[1,1], bbox(NY8)[2,2], labels="a)", cex=0.8)
plot(NY8, border="grey60", axes=TRUE)
points(TCE, pch=1, cex=0.7)
points(TCE, pch=3, cex=0.7)
text(coordinates(TCE), labels=as.character(TCE$name), cex=0.7,
     font=1, pos=c(4,1,4,1,4,4,4,2,3,4,2), offset=0.3)

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
nysar <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                  data = NY8, listw = NYlistw)
summary(nysar)
# conclude that (lambda, i.e., rho in our book) is different
# from zero according to a LRT
# exposure to TCE still not significant, but close enough
# that it shouldn't be discarded
# other two variables associated with higher incidence rates

# plot trend and stochastic components
# manually create color palettes for plot
NY8$sar_trend <- nysar$fit$signal_trend
NY8$sar_stochastic <- nysar$fit$signal_stochastic
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
tr_at <- seq(-1, 1.3, length.out=21)
tr_rds <- rds(sum(tr_at >= 0)*2)[-(1:(sum(tr_at >= 0)-sum(tr_at < 0)))]
st_at <- seq(-0.16, 0.39, length.out=21)
st_rds <- rds(sum(st_at >= 0)*2)[-(1:(sum(st_at >= 0)-sum(st_at < 0)))]
# create trellis graphic object for trend
tr_pl <- spplot(NY8, c("sar_trend"), at=tr_at, col="transparent",
                col.regions=tr_rds, main=list(label="Trend", cex=0.8))
# create trellis graphic object for stochastic component
st_pl <- spplot(NY8, c("sar_stochastic"), at=st_at, col="transparent",
                col.regions=st_rds, main=list(label="Stochastic", cex=0.8))
# plot the trellis objects
# the split arguments say, plot (the first plot) in position (1, 1) of a
# 2x1 matrix
# more is a logical specifying whether more plots will be following
plot(tr_pl, split = c(1,1,2,1), more = TRUE)
plot(st_pl, split = c(2,1,2,1), more = FALSE)

# changes weights from ls regression.  Weights are inversely proportional
# to population sizes in each region.
nylmw <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data = NY8,
            weights = POP8)
summary(nylmw)
# TCE exposure significant with expected sign!

# add residual from wls fit to NY8
NY8$lmwresid <- residuals(nylmw)

# create cool color palette (manually)
gry <- c(rev(brewer.pal(6, "Reds")[1:4]), colorRampPalette(brewer.pal(5, "Blues"))(9))
spplot(NY8, c("lmresid", "lmwresid"), col.regions=gry, col="transparent", lwd=0.5, at=seq(-2, 4.5, 0.5))

# check moran's after accounting for heterogeneous population sizes
lm.morantest(nylmw, NYlistw)
# model mispecification (rho != 0) in first test caused by
# heteroskedasticity more than to spatial autocorrelation?

# SAR model with different weights, that account for heterogeneous
# population sizes.  Fit sar model with different weights for Vv.  Specifically,
# Vv propto 1/ni for each region i
nysarw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                   data = NY8, listw = NYlistw, weights = POP8)
summary(nysarw)
# no traces of spatial autocorrelation after adjusting for the
# heterogeneous population sizes

# separate trend from stochastic component
NY8$sarw_trend <- nysarw$fit$signal_trend
NY8$sarw_stochastic <- nysarw$fit$signal_stochastic

# manually create color palettes for plot
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
tr_at <- seq(-1, 1.3, length.out=21)
tr_rds <- rds(sum(tr_at >= 0)*2)[-(1:(sum(tr_at >= 0)-sum(tr_at < 0)))]
st_at <- seq(-0.16, 0.39, length.out=21)
st_rds <- rds(sum(st_at >= 0)*2)[-(1:(sum(st_at >= 0)-sum(st_at < 0)))]

# create trellis graphic object for trend
tr_plw <- spplot(NY8, c("sarw_trend"), at=tr_at, col="transparent",
                col.regions=tr_rds, main=list(label="Trend", cex=0.8))
# create trellis graphic object for stochastic component
st_plw <- spplot(NY8, c("sarw_stochastic"), at=st_at, col="transparent",
                col.regions=st_rds, main=list(label="Stochastic", cex=0.8))
# plot the trellis objects
# the split arguments say, plot (the first plot) in position (1, 1) of a
# 2x1 matrix
# more is a logical specifying whether more plots will be following
plot(tr_plw, split = c(1,1,2,1), more = TRUE)
plot(st_plw, split = c(2,1,2,1), more = FALSE)
# resiual autocorrelation vanishes

# 2x1 matrix of plot of unweighted and weighted
# predicted spatial residuals for SAR models
plot(st_pl, split = c(1,1,2,1), more = TRUE)
plot(st_plw, split = c(2,1,2,1), more = FALSE)

# fit one-parameter car model
nycar <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                  data = NY8, family = "CAR", listw = NYlistw)
summary(nycar)
# the rho parameter is significant.  Suggests spatial
# autocorrelation

# fit car model with different weights for Vc.  Specifically,
# Vc propto 1/ni for each region i
nycarw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
                   data = NY8, family = "CAR", listw = NYlistw, weights = POP8)
summary(nycarw)
# spatial autocorrelation vanishes again

# extract fitted values from sarw and carw models
NY8$fitsarw = nysarw$fit$fitted.values
NY8$fitcarw = nycarw$fit$fitted.values

# plot results from sarw and carw (using two color palettes)
spplot(NY8, c("fitsarw", "fitcarw"), col.regions = fields::tim.colors(16))
library(colorRamps)
spplot(NY8, c("fitsarw", "fitcarw"), col.regions = blue2red(16))

# compare fitted values, manually
NY8$fitsarw - NY8$fitcarw

# overall conclusion
# we probably didn't need to account for spatial autocorrelation, assuming
# we accounted for varying population sizes.  Spatial patterns driven by
# population sizes, not true residual spatial autocorrelation.
