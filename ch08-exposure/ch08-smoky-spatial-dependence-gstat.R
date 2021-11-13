# install.packages(c("gstat", "geoR", "rgdal", "sp", "gridExtra", "sf", "colorspace"))
library(gstat) # for most of the work
library(sp)
library(rgdal)

### some notes
# a dataframe (df) can be turned into a SpatialPointsDataFrame
# using the command coordinates(df) <- c("x", "y") where
# x and y are the names of the x and y coordinate in df.

# the spplot function takes a SpatialPointsDataFrame or a
# the variogram function requires a gstat object, which
# can be created using the gstat function

# the predict function (to perform kriging) takes
# a gstat object that includes a variogram model

# the variogram functions for the gstat package parameterize
# the direction in degrees.  the variogram functions for the
# geoR package function the variogram functions in radians,
# where radians = degrees/180*pi

# projections for data and prediction locations must match

# turn smoky dataframe into SpatialPointsDataFrame by adding coordinates
load("./data/smoky.rda")
coordinates(smoky) <- c("longitude", "latitude")

# read polygon of data
poly = rgdal::readOGR("./data/smoky/smokypoly.shp")
proj4string(poly)
proj4string(smoky) #coordinate reference systems don't match!
proj4string(poly) = CRS(proj4string(smoky))

### Fig 8.4 create bubble plot of smoky pH
# place legend on right, change default colors with col.regions
spplot(smoky, "ph", key.space = "right", cuts = 10,
       col.regions = hcl.colors(11))

# create gstat object for further analysis
# formula ph ~ 1 assumes a constant mean over the spatial domain
gsmoky = gstat(id = "ph", formula = ph ~ 1, data = smoky)

### Fig 8.5
# construct variogram up to distance of 90 with bin widths of 4.5
variog1 = variogram(gsmoky, cutoff = 90, width = 4.5)
plot(variog1)

### Fig 8.6
# construct variogram up to distance of 120 with bin widths of 10
variog2 = variogram(gsmoky, cutoff = 120, width = 10)
plot(variog2)

# number of pairs used to estimate semivariance in each interval
variog2$np # gstat

### Table 8.1 fit variogram model using WRSS
# estimated exponential model with starting values c = .25, a = 30, c0 = .05
fitexp = fit.variogram(variog2, vgm(.25, "Exp", 30, .05),
                       fit.method = 2)
fitexp #c = 0.214, a = 22.55, c0 = 0.0138
# plot variogram with estimated exponential model
plot(variog2, fitexp, main = "WRSS exponential fit")
# wrss of fit
attr(fitexp, "SSErr")

### Table 8.2 # fit variogram model using WRSS
# estimated spherical model with starting values c = .25, a = 30, c0 = .05
fitsph = fit.variogram(variog2, vgm(.25, "Sph", 30, .05), fit.method = 2)
fitsph #c = 0.17, a = 63.33, c0 = 0.05
# plot variogram with estimated spherical model
plot(variog2, fitsph, main = "WRSS spherical fit")
# wrss of fit
attr(fitsph, "SSErr")
# the spherical model fits better since its WRSS = 15.9
# vs 19.9 for the exponential model

### fit variogram model using WRSS
# estimated matern model with starting values c = .25, a = 30, c0 = .05, smoothness = 1
fitmat = fit.variogram(variog2, vgm(.25, "Sph", 30, .05, kappa = 1), fit.method = 2,
                       fit.kappa = TRUE)
fitmat #c = 0.17, a = 63.33, c0 = 0.05
# plot variogram with estimated spherical model
plot(variog2, fitmat, main = "WRSS spherical fit")
# wrss of fit
attr(fitmat, "SSErr")
# the WRSS = 15.9

### Fig 8.11
# plot directional variogram
# direction are 70, 115, 160, 205
variog3 = variogram(gsmoky, alpha = c(70, 115, 160, 205))
plot(variog3)

### Fig 8.12
# show fit of anisotropic model
# the numbers come from output of the likfit function from the geoR package
# note that likfit parameterized a parameter in terms of aminor, not amajor
# so we have to convert it to amaj (amin * psiR).  We need to convert the ratio in anis
# to the minor/major ratio, i.e., 1/psiR
vgmaniso = vgm(.188, "Exp", 10.82*1.91, .0016, anis = c(70, 1/1.91))
plot(variog3, vgmaniso)