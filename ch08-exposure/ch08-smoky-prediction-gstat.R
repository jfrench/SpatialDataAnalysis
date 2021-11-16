# install.packages(c("gstat", "geoR", "rgdal", "sp", "gridExtra", "sf", "colorspace"))
library(gstat) # for most of the work
library(sp)
library(rgdal)
library(sf)
library(colorspace)

### some notes
# a dataframe (df) can be turned into a SpatialPointsDataFrame
# using the command coordinates(df) <- c("x", "y") where
# x and y are the names of the x and y coordinate in df.

# the predict function (to perform kriging) takes
# a gstat object that includes a variogram model

# the variogram functions for the gstat package parameterize
# the direction in degrees.  the variogram functions for the

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

# make variograms from book
# anisotropic, then filtered anisotropic
# c = 0.2725, model = exponential, amaj = 36.25, direction of amajor is 70 degrees,
# ratio of amin/amaj = 16.93/36.25
# C0 = 0.0325, the a chould be amaj, and the ratio in anis should be amin/amaj
# v1 is the variogram model for standard kriging
# v2 is the variogram model for filtered kriging
v1 = vgm(0.2725, "Exp", 36.25, nugget = 0.0325, anis = c(70, 16.93/36.25))
v2 = vgm(0.2725, "Exp", 36.25, Err = 0.0325, anis = c(70, 16.93/36.25))

### Fig 8.14
# create prediction grid
grid = spsample(poly, n = 1600, type = "regular") # grid of points within polygon
coordnames(grid) = c("longitude", "latitude") # coordinate names have to match original data
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
       col.regions = hcl.colors(64), cuts = 63,
       main = "ok predictions")


# plot kriging variance from ok model
spplot(ok, "ph.var", colorkey = TRUE,
       col.regions = hcl.colors(64), cuts = 63,
       main = "ok kriging variance")
# plot predictions from fok model
spplot(fok, "ph.pred", colorkey = TRUE,
       col.regions = hcl.colors(64), cuts = 63, main = "filtered ok predictions")
# plot kriging variance from fok model
spplot(fok, "ph.var", colorkey = TRUE,
       col.regions = hcl.colors(64), cuts = 63, main = "filtered ok kriging variance")
# plot predictions from idw
spplot(idp, "var1.pred", colorkey = TRUE,
       col.regions = hcl.colors(64), cuts = 63, main = "idw predictions")

# using ggplot2 (perhaps not the most efficient way)
library(sf)
sfok <- sf::st_as_sf(ok) # keep points as columns
# add u- and v-positions to sfok data frame
sfok = cbind(sfok, sf::st_coordinates(sfok))
# get u- and v-position of original data
coordsdf = as.data.frame(coordinates(smoky))
library(ggplot2)
# plot predictions
ggplot() + # start ggplot
        geom_tile(aes(x = X, y = Y, fill = ph.pred), data = sfok) + # create heat map
        scale_fill_viridis_c() +
        coord_fixed() +  # change coordinate reference system
        geom_point(aes(x = longitude, y = latitude), data = coordsdf) # add observed points

# plot kriging variance
ggplot() + # start ggplot
  geom_tile(aes(x = X, y = Y, fill = ph.var), data = sfok) + # create heat map
  scale_fill_viridis_c() +
  coord_fixed() +  # change coordinate reference system
  geom_point(aes(x = longitude, y = latitude), data = coordsdf) # add observed points

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
                 col.regions = hcl.colors(64), at = cut,
                 main = "probability map of Pr(ph > 8) ordinary")
ifkplot = spplot(ifk, "ph.pred", colorkey = TRUE,
                 col.regions = hcl.colors(64), at = cut,
                 main = "probability map of Pr(ph > 8) filtered")
library(gridExtra)
grid.arrange(iokplot, ifkplot, ncol = 2)

#### example of conditional simulation
# grid of points, fairly coarse.
# don't converted to gridded data
grid2 = spsample(poly, n = 100, type = "regular")
# don't convert grid2 to a grid
ok_sim = predict(ganiso1, newdata = grid2, nsim = 100)
spplot(ok_sim, "sim1", col.regions = hcl.colors(11), cuts = 10)

gridded(grid2) = TRUE
ok_sim = predict(ganiso1, newdata = grid2, nsim = 100)
spplot(ok_sim, "sim1", col.regions = hcl.colors(11))

coordnames(grid2) = c("longitude", "latitude") # coordinate names have to match original data
gridded(grid2) = TRUE # turn into grid for better plotting!
ok_sim = predict(ganiso1, newdata = grid2, nsim = 100)
# plot some of the surfaced.  Colors not on same scale
spplot(ok_sim, "sim1", col.regions = hcl.colors(11), cuts = 10)
spplot(ok_sim, "sim2", col.regions = hcl.colors(11), cuts = 10)
spplot(ok_sim, "sim3", col.regions = hcl.colors(11), cuts = 10)

range_sim = range(data.frame(attr(ok_sim, "data"))) # get range of all conditional simulations
cut = seq(min(range_sim), max(range_sim), len = 63) # for consistent coloring of graphics
# construct plots with consistent coloring
sim1plot = spplot(ok_sim, "sim1",
                  col.regions = hcl.colors(64), at = cut,
                  main = "heat map of 1st simulation")
sim2plot = spplot(ok_sim, "sim2",
                  col.regions = hcl.colors(64), at = cut,
                  main = "heat map of 2nd simulation")
sim3plot = spplot(ok_sim, "sim3",
                  col.regions = hcl.colors(64), at = cut,
                  main = "heat map of 3rd simulation")
grid.arrange(sim1plot, sim2plot, sim3plot, ncol = 3)

## using ggplot2
sf_sim <- sf::st_as_sf(ok_sim) # keep points as columns
# add u- and v-positions to sfok data frame
sf_sim = cbind(sf_sim, sf::st_coordinates(sf_sim))
sf_sim <- sf_sim[, c("X", "Y", "sim1", "sim2", "sim3")]

# need to make data long (all z values are in 1 columns)
sf_sim_long <- sf_sim |> tidyr::pivot_longer(cols = c("sim1", "sim2", "sim3"))
ggplot(data = sf_sim_long) + geom_tile(aes(x = X, y = Y, fill = value)) +
  facet_grid(~ name) + scale_fill_viridis_c()



