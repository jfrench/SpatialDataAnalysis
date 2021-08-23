# install.packages(c("gstat", "geoR", "fields", "maptools", "gridExtra"))
library(gstat) # for most of the work
library(geoR) # for estimating anisotropy
library(sp)

# set working directory
setwd("~/OneDrive - The University of Colorado Denver/Teaching/Math6384/Data")
### some notes
# a dataframe (df) can be turned into a SpatialPointsDataFrame
# using the command coordinates(df) <- c("x", "y") where
# x and y are the names of the x and y coordinate in df.

# the spplot function takes a SpatialPointsDataFrame or a
# the variogram function requires a gstat object, which
# can be created using the gstat function

# the predict function (to perform kriging) takes
# a gstat object that includes a variogram model

# the variog and likfit functions (from the geoR package)
# required geodata objects, which are created using the
# as.geodata function.

# the variogram functions for the gstat package parameterize
# the direction in degrees.  the variogram functions for the
# geoR package function the variogram functions in radians,
# where radians = degrees/180*pi

# projections for data and prediction locations must match

# turn smoky dataframe into SpatialPointsDataFrame by adding coordinates
load("smoky.rda")
coordinates(smoky) <- c("easting", "northing")

# read polygon of data
poly = rgdal::readOGR("smoky/smokypoly.shp")
proj4string(poly)
proj4string(smoky) #coordinate reference systems don't match!
proj4string(poly) = CRS(proj4string(smoky))

### Fig 8.4 create bubble plot of smoky pH
# place legend on right, change default colors with col.regions
spplot(smoky, "ph", key.space = "right", cuts = 10,
       col.regions = colorspace::terrain_hcl(11))

# create gstat object for further analysis
# formula ph ~ 1 assumes a constant mean over the spatial domain
gsmoky = gstat(id = "ph", formula = ph ~ 1, data = smoky)

# create geodata object for geoR
geosmoky = as.geodata(cbind(smoky$easting, smoky$northing, smoky$ph))

### Fig 8.5
# construct variogram up to distance of 90 with bin widths of 4.5
variog1 = variogram(gsmoky, cutoff = 90, width = 4.5)
plot(variog1)

# same plot using geoR
variog1b = variog(geosmoky, max.dist = 90, uvec = 20)
plot(variog1b)

### Fig 8.6
# construct variogram up to distance of 120 with bin widths of 10
variog2 = variogram(gsmoky, cutoff = 120, width = 10)
plot(variog2)

# same plot using geoR
variog2b = variog(geosmoky, max.dist = 120, uvec = 12)
plot(variog2b)

# number of pairs used to estimate semivariance in each interval
variog2$np # gstat
variog2b$n # geoR

### Table 8.1 fit variogram model using WRSS
# estimated exponential model with starting values c = .25, a = 30, c0 = .05
fitexp = fit.variogram(variog2, vgm(.25, "Exp", 30, .05),
                       fit.method = 2)
fitexp #c = 0.214, a = 22.55, c0 = 0.0138
# plot variogram with estimated exponential model
plot(variog2, fitexp, main = "WRSS exponential fit")

# same thing using geoR (numbers a little different because variogram a little different)
fitexpb = variofit(variog2b, ini.cov.pars = c(.25, 30),
                   nugget = 0.05, cov.model = "exponential",
                   weights = "cressie")
fitexpb
plot(variog2b)
lines(fitexpb)

### Table 8.2 # fit variogram model using WRSS
# estimated spherical model with starting values c = .25, a = 30, c0 = .05
fitsph = fit.variogram(variog2, vgm(.25, "Sph", 30, .05), fit.method = 2)
fitsph #c = 0.17, a = 63.33, c0 = 0.05
# plot variogram with estimated spherical model
plot(variog2, fitsph, main = "WRSS spherical fit")

# same thing using geoR (numbers different because variograms a little different)
fitsphb = variofit(variog2b, ini.cov.pars = c(.25, 30), nugget = 0.05, cov.model = "spherical", weights = "cressie")
fitsphb
plot(variog2b)
lines(fitsphb)

# fit variogram model using Restricted Maximum Likelihood
# gstat isn't very good at this, so we'll use the likfit function
# from the geoR package

# use REML to estimate the parameters of geosmoky for an exponential model
# with starting values c = .25, a = 30, c0 = 0.05
lfit = likfit(geosmoky, ini.cov.pars = c(.25, 30),
              nugget = 0.05,
              cov.model = "exponential",
              lik.method = "REML")
lfit # c = .1929, a = 10.97, c0 = 0
# plot estimated model
plot(variog2, vgm(.1929, "Exp", 10.97, 0), main = "REML exponential fit")

### Fig 8.11
# plot directional variogram
# direction are 70, 115, 160, 205
variog3 = variogram(gsmoky, alpha = c(70, 115, 160, 205))
plot(variog3)

# same (general) thing using geoR.  Angles must be between 0 and 180.  I get 25 from 205 - 180
# geoR specifies direction in radians, which is degrees/180*pi
variog3b = variog4(geosmoky, direction = c(25, 70, 115, 160)/180*pi, max.dist = 120)
plot(variog3b)

# fit REML model for anisotropy
# psiA controls the anisotropy angle
# psiR controls the ratio of the large a (amaj) to the smallest a (amin)
# note that psiA must be converted from degrees to radians using the formula radians = degrees/180*pi
# by default psiA and psiR are considered FIXED (not estimates)
# in model 1, both are fixed, in model 2, the psiA is fixed, and in model 3, neither are fixed
lfit0 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05, cov.model = "exponential", lik.method = "REML") # omnidirectional
lfit1 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05, cov.model = "exponential", lik.method = "REML", psiA = 70/180*pi, psiR = 3)
lfit2 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05, cov.model = "exponential", lik.method = "REML", psiA = 70/180*pi, psiR = 3, fix.psiR = FALSE)
lfit3 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05, cov.model = "exponential", lik.method = "REML", psiA = 70/180*pi, psiR = 3, fix.psiR = FALSE, fix.psiA = FALSE)
# fit matern variogram model.  kappa is the smoothness parameter.  we allow it to be estimated.
lfit4 = likfit(geosmoky, ini.cov.pars = c(.25, 30), nugget = 0.05, cov.model = "matern", lik.method = "REML", kappa = 1, fix.kappa = FALSE)


# compare AIC/BIC values.  Model 3 is best in terms of AIC
# but we'll use model 2
summary(lfit1)
summary(lfit2) # c = .188, amin = 10.82, c0= 0.0016, amaj/amin = 1.91, degrees = 70
summary(lfit3)
summary(lfit4)

### Fig 8.12
# show fit of anisotropic model
# note that likfit parameterized a parameter in terms of aminor, not amajor
# so we have to convert it to amaj (amin * psiR).  We need to convert the ratio in anis
# to the minor/major ratio, i.e., 1/psiR
vgmaniso = vgm(.188, "Exp", 10.82*1.91, .0016, anis = c(70, 1/1.91))
plot(variog3, vgmaniso)

# variogram for filtered data

# make variograms from book
# anisotropic, then filtered anisotropic
# the a chould be amaj, and the ratio in anis should be amin/amaj
v1 = vgm(0.2725, "Exp", 36.25, nugget = 0.0325, anis = c(70, 16.93/36.25))
v2 = vgm(0.2725, "Exp", 36.25, Err = 0.0325, anis = c(70, 16.93/36.25))

### Fig 8.14
# create prediction grid
grid = spsample(poly, n = 1600, type = "regular") # grid of points within polygon
coordnames(grid) = c("easting", "northing") # coordinate names have to match original data
gridded(grid) = TRUE # turn into grid for better plotting!

# gstat objects with anisotropic variogram models for kriging
ganiso1 = gstat(id = "ph", formula = ph ~ 1, data = smoky, model = v1)
ganiso2 = gstat(id = "ph", formula = ph ~ 1, data = smoky, model = v2)

# kriging requires a gstat object, as create above
# idw only needs a data frame (or SpatialPointsDataFrame)
# ordinary kriging
ok = predict(ganiso1, newdata = grid)
# filtered ordinary kriging
fok = predict(ganiso2, newdata = grid)
# inverse distance weighting prediction
idp = idw(ph ~ 1, smoky, newdata = grid)

# plot predictions from ok model
spplot(ok, "ph.pred", colorkey = TRUE,
       col.regions = colorspace::terrain_hcl(64), cuts = 63,
       main = "ok predictions")


# plot kriging variance from ok model
spplot(ok, "ph.var", colorkey = TRUE,
       col.regions = colorspace::terrain_hcl(64), cuts = 63,
       main = "ok kriging variance")
# plot predictions from fok model
spplot(fok, "ph.pred", colorkey = TRUE,
       col.regions = colorspace::terrain_hcl(64), cuts = 63, main = "filtered ok predictions")
# plot kriging variance from fok model
spplot(fok, "ph.var", colorkey = TRUE,
       col.regions = colorspace::terrain_hcl(64), cuts = 63, main = "filtered ok kriging variance")
# plot predictions from idw
spplot(idp, "var1.pred", colorkey = TRUE,
       col.regions = colorspace::terrain_hcl(64), cuts = 63, main = "idw predictions")

# using ggplot2 (perhaps not the most efficient way)
library(sf)
# add u- and v-positions to sfok data frame
sfok2 = cbind(sfok, sf::st_coordinates(sfok))
# get u- and v-position of original data
coordsdf = as.data.frame(coordinates(smoky))
library(colorspace)
library(ggplot2)
# plot predictions
ggplot() + # start ggplot
        geom_tile(aes(x = X, y = Y, fill = ph.pred), data = sfok2) + # create heat map
        scale_fill_gradientn(colours = colorspace::terrain_hcl(64)) + # change color scale
        coord_fixed() +  # change coordinate reference system
        geom_point(aes(x = easting, y = northing), data = coordsdf) # add observed points

# plot kriging variance
ggplot() + # start ggplot
        geom_tile(aes(x = X, y = Y, fill = ph.var), data = sfok2) + # create heat map
        scale_fill_gradientn(colours = colorspace::terrain_hcl(64)) + # change color scale
        coord_fixed() +  # change coordinate reference system
        geom_point(aes(x = easting, y = northing), data = coordsdf) # add observed points

#### Example of indicator kriging
# indicator kriging for ph > 8
# ordinary indicator kriging
ganiso_i1 = gstat(id = "ph", formula = (ph > 8 + 0) ~ 1, data = smoky, model = v1)
iok = predict(ganiso_i1, newdata = grid)
# filtered indicator kriging
ganiso_i2 = gstat(id = "ph", formula = (ph > 8 + 0) ~ 1, data = smoky, model = v2)
ifk = predict(ganiso_i2, newdata = grid)

# compare maps of exceedance probabilities (make sure on same scale)
range(iok$ph.pred) # negative probabilities!
allp = c(iok$ph.pred, ifk$ph.pred)
cut = seq(min(allp), 1, len = 63) # for consistent coloring of graphics
iokplot = spplot(iok, "ph.pred", colorkey = TRUE,
                 col.regions = colorspace::terrain_hcl(64), at = cut,
                 main = "probability map of Pr(ph > 8) ordinary")
ifkplot = spplot(ifk, "ph.pred", colorkey = TRUE,
                 col.regions = colorspace::terrain_hcl(64), at = cut,
                 main = "probability map of Pr(ph > 8) filtered")
library(gridExtra)
grid.arrange(iokplot, ifkplot, ncol = 2)

#### example of conditional simulation
# grid of points, fairly coarse.
# don't converted to gridded data
grid2 = spsample(poly, n = 100, type = "regular")
# don't convert grid2 to a grid
ok_sim = predict(ganiso1, newdata = grid2, nsim = 100)
spplot(ok_sim, "sim1")

gridded(grid2) = TRUE
ok_sim = predict(ganiso1, newdata = grid2, nsim = 100)
spplot(ok_sim, "sim1")

coordnames(grid2) = c("easting", "northing") # coordinate names have to match original data
gridded(grid2) = TRUE # turn into grid for better plotting!
ok_sim = predict(ganiso1, newdata = grid2, nsim = 100)
# plot some of the surfaced.  Colors not on same scale
spplot(ok_sim, "sim1")
spplot(ok_sim, "sim2")
spplot(ok_sim, "sim3")

range_sim = range(data.frame(attr(ok_sim, "data"))) # get range of all conditional simulations
cut = seq(min(range_sim), max(range_sim), len = 63) # for consistent coloring of graphics
# construct plots with consistent coloring
sim1plot = spplot(ok_sim, "sim1",
                  col.regions = colorspace::terrain_hcl(64), at = cut,
                  main = "heat map of 1st simulation")
sim2plot = spplot(ok_sim, "sim2",
                  col.regions = colorspace::terrain_hcl(64), at = cut,
                  main = "heat map of 2nd simulation")
sim3plot = spplot(ok_sim, "sim3",
                  col.regions = colorspace::terrain_hcl(64), at = cut,
                  main = "heat map of 3rd simulation")
grid.arrange(sim1plot, sim2plot, sim3plot, ncol = 3)




