library(gstat) # for most of the work
library(sf)

load("./data/smoky.rda")

# turn smoky data frame into sf data frame
# remove = TRUE means keep coordinate columns as
# variables in data frame
smoky <- sf::st_as_sf(smoky,
                      coords = c("longitude", "latitude"),
                      remove = FALSE)

plot(smoky["ph"], pch = 20, pal = hcl.colors)

# read polygon of data
poly = sf::st_read("./data/smoky/smokypoly.shp")
sf::st_crs(poly)
sf::st_crs(smoky) #coordinate reference systems don't match!
sf::st_crs(poly) = sf::st_crs(smoky)

# verify that they are in fact the same crs
plot(st_geometry(poly))
points(st_coordinates(smoky))

# create prediction grid
pgrid = st_sample(poly, size = 1600, type = "regular")
# extract coordinates from pgrid
pcoords = st_coordinates(pgrid)
# need to add longitude and latitude as variables in sf data frame
pdf = data.frame(longitude = pcoords[,1],
                 latitude = pcoords[,2])
pdf = st_as_sf(pdf, coords = c("longitude", "latitude"),
               remove = FALSE)

# universal kriging in gstat
# notice that formula includes predictor variables
uksmoky = gstat(id = "ph",
                formula = ph ~ longitude + latitude,
                data = smoky)
# create directional variogram for residuals to see if we have anisotropy
variog4 = variogram(uksmoky, cutoff = 90,
                    alpha = c(70, 115, 160, 205))
plot(variog4) # no evidence of anisotropy

# create omnidirectional variogram for uk
variog = variogram(uksmoky, cutoff = 90)
v = fit.variogram(variog, vgm(.3, "Exp", 20, 0), fit.method = 2)
v # see estimates
plot(variog, v) # fits well

# add variogram model to uksmoky
uksmoky = gstat(id = "ph",
                formula = ph ~ longitude + latitude,
                data = smoky,
                model = v)

# make universal kriging predictions on grid
uk = predict(uksmoky, pdf)

# plot prediction from uk
plot(uk["ph.pred"], pal = hcl.colors)
plot(uk["ph.var"], pal = hcl.colors)

# evaluate covariance matrix for observed data
# determine distances
d = as.matrix(dist(st_coordinates(smoky)))
# determine estimated covariance matrix
C = .131 * exp(-d/18.97) + 0.043 * diag(nrow(st_coordinates(smoky)))
# create X matrix
X = cbind(1, st_coordinates(smoky))
# determine betahat gls
betahat_gls = solve(crossprod(X, solve(C, X)), t(X) %*% solve(C, smoky$ph))
# determine residuals
r = smoky$ph - X %*% betahat_gls
# add residuals to smoky df
smoky$r = r
# create gstat object for smoky residuals
rsmoky = gstat(id = "r",
               formula = r ~ 1,
               data = smoky,
               model = v)
# make ordinary kriging predictions of residual on grid
rok = predict(rsmoky, pdf)
# add trend back into ordinary kriging predictions
rhat = cbind(1, st_coordinates(pdf)) %*% betahat_gls + rok$r.pred

# predictions the same (subject to rounding)
range(rhat - uk$ph.pred)
head(cbind(rhat, uk$ph.pred))

# add ok variance to uk object for comparison
uk$ph.var.ok <- rok$r.var

# plot kriging variances
plot(uk[c("ph.var", "ph.var.ok")],
     pch = 20,
     pal = hcl.colors, key.pos = -1)

# difference in variances
uk$ph.var.diff <- uk$ph.var - uk$ph.var.ok

plot(uk["ph.var.diff"], pal = hcl.colors, pch = 20)