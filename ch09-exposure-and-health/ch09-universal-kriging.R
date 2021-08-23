library(gstat) # for most of the work
library(sp) # for plotting

setwd("~/OneDrive - The University of Colorado Denver/Teaching/Math6384/Data")
load("smoky.rda")

# turn smoky dataframe into SpatialPointsDataFrame by adding coordinates
coordinates(smoky) <- c("easting", "northing")

### create bubble plot of smoky pH
# place legend on right, change default colors with col.regions
spplot(smoky, "ph", key.space = "right", cuts = 10, col.regions = topo.colors(11))

# create prediction grid
library(maptools) # to read shape polygon
# read polygon of data
poly = readShapePoly("smoky/smokypoly.shp")

grid = spsample(poly, n = 1600, type = "regular") # grid of points within polygon
coordnames(grid) = c("easting", "northing") # coordinate names have to match original data
gridded(grid) = TRUE # turn into grid for better plotting!

# universal kriging in gstat
# notice that formula includes predictor variables
uksmoky = gstat(id = "ph", formula = ph ~ easting + northing, data = smoky)
# create directional variogram for residuals to see if we have anisotropy
variog4 = variogram(uksmoky, cutoff = 90, alpha = c(70, 115, 160, 205))
plot(variog4) # no evidence of anisotropy

# create omnidirectional variogram for uk
variog = variogram(uksmoky, cutoff = 90)
v = fit.variogram(variog, vgm(.3, "Exp", 20, 0), fit.method = 2)
v # see estimates
plot(variog, v) # fits well

# add variogram model to uksmoky
uksmoky = gstat(id = "ph", formula = ph ~ easting + northing, data = smoky, model = v)

# make universal kriging predictions on grid
uk = predict(uksmoky, grid)

# plot prediction from uk
spplot(uk, "ph.pred", colorkey = TRUE, col.regions = fields::tim.colors(64), cuts = 63, main = "uk predictions")
# plot kriging variance from uk
spplot(uk, "ph.var", colorkey = TRUE, col.regions = fields::tim.colors(64), cuts = 63)

# evaluate covariance matrix for observed data
# determine distances
d = as.matrix(dist(coordinates(smoky)))
# determine estimated covariance matrix
C = .131 * exp(-d/18.97) + 0.043 * diag(nrow(coordinates(smoky)))
# create X matrix
X = cbind(1, coordinates(smoky))
# determine betahat gls
betahat_gls = solve(crossprod(X, solve(C, X)), t(X) %*% solve(C, smoky$ph))
# determine residuals
r = smoky$ph - X %*% betahat_gls
# add residuals to smoky df
smoky$r = r
# create gstat object for smoky residuals
rsmoky = gstat(id = "r", formula = r ~ 1, data = smoky, model = v)
# make ordinary kriging predictions of residual on grid
rok = predict(rsmoky, grid)
# add trend back into ordinary kriging predictions
rhat = cbind(1, coordinates(grid)) %*% betahat_gls + rok$r.pred

# predictions the same (subject to rounding)
range(rhat - uk$ph.pred)
head(cbind(rhat, uk$ph.pred))

# plot ordinary and universal kriging variances on one plot
cut = seq(0, max(uk$ph.var), len = 63) # for consistent coloring of graphics
# construct plots with consistent coloring
okvar = spplot(rok, "r.var", col.regions = colorspace::terrain_hcl(64),
               at = cut, main = "ordinary kriging variance")
ukvar = spplot(uk, "ph.var", col.regions = colorspace::terrain_hcl(64),
               at = cut, main = "universal kriging variance")
library(gridExtra)
grid.arrange(okvar, ukvar, ncol = 2)
# compare kriging variances
head(cbind(rok$r.var, uk$ph.var))
# determine proportions of kriging variances that are
# smaller for ok than uk
sum(rok$r.var <= uk$ph.var)/length(rok$r.var)
