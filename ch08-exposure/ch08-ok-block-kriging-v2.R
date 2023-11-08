# install.packages(c("gstat", "geoR", "fields", "maptools", "gridExtra"))
library(gstat) # for most of the work
library(sf)
library(ggplot2)

# load colorado data
data(co, package = "gear")
co <- st_as_sf(co, coords = c("easting", "northing"))
plot(co["Al"], pal = hcl.colors, pch = 20)

# read colorado zctas
# set working directory
zctas = st_read("./data/co_zcta/Colorado_ZIP_Code_Tabulation_Areas_ZCTA.shp")
plot(zctas, axes = TRUE)
# convert longitude/latitude to UTM coordinates zone 13
# this is approximately EPSG 26913
# https://epsg.io/26913
zctas = st_transform(zctas, 26913)
plot(zctas["GEOID10"], axes = TRUE)

# ensure coordinate reference systems are the same
st_crs(co) <- st_crs(zctas)

# ensure crs line up visually
plot(st_geometry(zctas), axes = TRUE)
points(st_coordinates(co), col = "blue", pch = 20)

# estimate omnidirectional semivariogram
# for aluminum (for simplicity, we won't
# examine anisotropy, but we should)
gco = gstat(id = "Al",
            formula = Al ~ 1,
            data = co)

# crude fit
vco = variogram(gco, cutoff = 350000)
fit.vco = fit.variogram(vco, vgm(psill = 0.5, model = "Sph",
                                 range = 3e5, nugget = 1),
                        fit.method = 2)
plot(vco, fit.vco)

# create gstat object with estimated semivariogram
gco = gstat(id = "Al",
            formula = Al ~ 1,
            data = co,
            model = fit.vco)

# point kriging prediction at centroids of zcta blocks
# get centroids -> turn into coordinates -> turn into data.frame
zcta_centroids <- data.frame(st_coordinates((st_centroid(zctas))))
# change column names
names(zcta_centroids) <- c("easting", "northing")
# convert to sf data frame
zcta_centroids <- st_as_sf(zcta_centroids,
                           coords = c("easting", "northing"))
# make sure sf data frames have the same crs
st_crs(zcta_centroids) <- st_crs(zctas)

# point kriging
p_point = predict(gco, newdata = zcta_centroids)

# block kriging
p_block = predict(gco, newdata = zctas)

# combine predictions into single data frame
p_point_block <- rbind(
  cbind(p_point[c("Al.pred", "Al.var")], type = "point"),
  cbind(p_block[c("Al.pred", "Al.var")], type = "block"))

# plot of predictions
ggplot(p_point_block) +
  geom_sf(aes(fill = Al.pred, color = Al.pred)) +
  facet_wrap(~ type) +
  scale_fill_viridis_c() +
  scale_color_viridis_c()

# plot of kriging variance
ggplot(p_point_block) +
  geom_sf(aes(fill = Al.var, color = Al.var)) +
  facet_wrap(~ type) +
  scale_fill_viridis_c() +
  scale_color_viridis_c()



