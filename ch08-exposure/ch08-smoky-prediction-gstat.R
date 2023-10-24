# install.packages(c("gstat", "geoR", "rgdal", "sp", "gridExtra", "sf", "colorspace"))
library(gstat) # for most of the work
library(sf)

### some notes
# a dataframe (df) can be turned into an sf
# data frame using st_as_sf

# the predict function (to perform kriging) takes
# a gstat object that includes a variogram model

# the variogram functions for the gstat package parameterize
# the direction in degrees.  the variogram functions for the

# projections for data and prediction locations must match

# turn smoky dataframe into sf data frame
load("./data/smoky.rda")
smoky <- sf::st_as_sf(smoky,
                      coords = c("longitude", "latitude"))

# read polygon of data
poly = sf::st_read("./data/smoky/smokypoly.shp")
sf::st_crs(poly)
sf::st_crs(smoky) #coordinate reference systems don't match!
sf::st_crs(poly) = sf::st_crs(smoky)


### Fig 8.4 create bubble plot of smoky pH
# change default colors with pal
# change point style with pch
plot(smoky["ph"], nbreaks = 10, pal = hcl.colors, pch = 20)

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
# Err means the nugget is measurement error, indicating we want filtered kriging
v1 = vgm(0.2725, "Exp", 36.25, nugget = 0.0325, anis = c(70, 16.93/36.25))
v2 = vgm(0.2725, "Exp", 36.25, Err = 0.0325, anis = c(70, 16.93/36.25))

### Fig 8.14
# create prediction grid
grid = st_sample(poly, size = 1600, type = "regular") # grid of points within polygon

# gstat objects with anisotropic variogram models for kriging
ganiso1 = gstat(id = "ph", formula = ph ~ 1, data = smoky,
                model = v1)
ganiso2 = gstat(id = "ph", formula = ph ~ 1, data = smoky,
                model = v2)

# kriging requires a gstat object
# idw only needs a data frame
# ordinary kriging
ok = predict(ganiso1, newdata = grid)
# filtered ordinary kriging
fok = predict(ganiso2, newdata = grid)
# inverse distance weighting prediction
idp = idw(ph ~ 1, smoky, newdata = grid)

# plot predictions from ok model
library(ggplot2)

# plot kriging predictions from ok model
ggplot(ok) +
  geom_sf(aes(col = ph.pred)) +
  scale_color_viridis_c()

# plot kriging variance from ok model
ggplot(ok) +
  geom_sf(aes(col = ph.var)) +
  scale_color_viridis_c()

# plot predictions from fok model
ggplot(fok) +
  geom_sf(aes(col = ph.pred)) +
  scale_color_viridis_c()

# plot kriging variance from fok model
ggplot(fok) +
  geom_sf(aes(col = ph.var)) +
  scale_color_viridis_c()

# plot predictions from idw
ggplot(idp) +
  geom_sf(aes(col = var1.pred)) +
  scale_color_viridis_c()

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

# predictions the same at unsampled locations!
all.equal(iok$ph.pred, ifk$ph.pred)

# create data frame to create plots side by side
iksf = rbind(cbind(iok, type = "unfiltered"),
             cbind(ifk, type = "filtered"))

# plot predictions from indicator kriging model
ggplot(iksf) + geom_sf(aes(col = ph.pred)) +
  scale_color_viridis_c() +
  facet_wrap(~ type) + ggtitle("indicator kriging predictions")

# plot kriging variance from indicator kriging model
ggplot(iksf) + geom_sf(aes(col = ph.var)) +
  scale_color_viridis_c() +
  facet_wrap(~ type) + ggtitle("indicator kriging variance")

#### example of conditional simulation
# grid of points, fairly coarse.
grid2 = st_sample(poly, size = 100, type = "regular")
# simulate at gridded locations
ok_sim = predict(ganiso1, newdata = grid2, nsim = 100)

# plot conditional simulations
ggplot(ok_sim) + geom_sf(aes(col = sim1)) +
  scale_color_viridis_c()

ggplot(ok_sim) + geom_sf(aes(col = sim2)) +
  scale_color_viridis_c()

ggplot(ok_sim) + geom_sf(aes(col = sim3)) +
  scale_color_viridis_c()

# plot conditional simulations in one plot

# need to make data long (all z values are in 1 column)
# also need x, y coordinates as variable
# the first part of the code creates a single column containing
# sim1, sim2, sim3
# the second part adds the associated X, Y coordinates as
# variables in the data frame
ok_sim_long <-
  ok_sim[c("sim1", "sim2", "sim3")] |>
  tidyr::pivot_longer(cols = c("sim1", "sim2", "sim3")) |>
  tibble::add_column(X = rep(st_coordinates(ok_sim)[,1], times = 3),
                     Y = rep(st_coordinates(ok_sim)[,2], times = 3))

# plot 3 simulations simultaneously using plot.sf
plot(ok_sim[c("sim1", "sim2", "sim3")], pal = hcl.colors,
     pch = 20, key.pos = 4)

# plot 3 simulations simultaneously using geom_sf
ggplot(ok_sim_long) + geom_sf(aes(col = value)) +
  facet_wrap(~name) + scale_color_viridis_c()

# plot rasterized version of 3 simulations using geom_tile
ggplot(ok_sim_long) +
  geom_raster(aes(x = X, y = Y, fill = value)) +
  facet_wrap(~name) + scale_fill_viridis_c()