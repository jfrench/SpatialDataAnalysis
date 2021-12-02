library(gear) # for most of the work
library(sf) # for plotting

load("./data/smoky.rda")
smoky_points <- st_as_sf(smoky, coords = c("longitude", "latitude"))

### create bubble plot of smoky pH
# place legend on right, change default colors with col.regions
plot(smoky_points["ph"], nbreaks = 11, pal = hcl.colors, pch = 19)

# read polygon of data
smoky_poly = st_read("./data/smoky/smokypoly.shp")
smoky_poly <- st_set_crs(smoky_poly, NA)
# ensure CRS matches
smoky_points <- st_set_crs(smoky_points, st_crs(smoky_poly))

# universal kriging in gear (using maximum likelihood)

# construct empirical semivariogram to estimate initial parameters
smoky_evgram <- evgram(ph ~ longitude + latitude, data = smoky,
                       coordnames = ~ longitude + latitude)
plot(smoky_evgram)
lattice::xyplot(smoky_evgram)
ggplot2::autoplot(smoky_evgram)

# fit geostatistical linear model to data with estimates
smoky_geolmod <- geolm(ph ~ longitude + latitude, data = smoky,
                       mod = cmod_std("exponential",
                                      psill = 0.15,
                                      r = 20,
                                      evar = 0.05),
                       coordnames = ~ longitude + latitude)
# update geostatistical linear model with REML estimates (isotropic model)
smoky_geolmod_iso <- estimate(smoky_geolmod, reml = TRUE)
smoky_geolmod_iso
# estimated covariance model parameters
smoky_geolmod_iso$mod

# update geostatistical linear model with REML estimates for a
# geometric anisotropic model estimating both the direction of minimum spatial
# dependence and the ratio of minimum spatial dependence to maximum spatial
# dependence
smoky_geolmod_aniso <- estimate(smoky_geolmod,
                                reml = TRUE,
                                est_angle = TRUE,
                                est_ratio = TRUE)
smoky_geolmod_aniso
# estimated covariance model parameters
smoky_geolmod_aniso$mod

# create prediction grid
pgrid <- st_sample(smoky_poly, 1000, type = "regular", exact = FALSE)
plot(pgrid)
# extract grid coordinates
ppoints <- st_coordinates(pgrid)
# convert grid to data frame with proper names
ppoints <- as.data.frame(ppoints)
colnames(ppoints) <- c("longitude", "latitude")

# predict at prediction locations (isotropic model)
piso <- predict(smoky_geolmod_iso, newdata = ppoints, return_type = "sf")
# bubble plot
plot(piso, pal = hcl.colors)

# predict at prediction locations (anisotropic model)
paniso <- predict(smoky_geolmod_aniso, newdata = ppoints, return_type = "sf")
# bubble plot
plot(paniso, pal = hcl.colors)
