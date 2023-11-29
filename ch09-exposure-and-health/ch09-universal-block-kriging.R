# install.packages(c("gstat", "geoR", "fields", "maptools", "gridExtra"))
library(gstat) # for most of the work
library(sf)
library(splines) # for ns
library(ggplot2) # for closing plots

# load colorado data
data(co, package = "gear")
co <- st_as_sf(co, coords = c("easting", "northing"),
               remove = FALSE)
plot(co["Al"], pal = hcl.colors, pch = 20)

# read colorado zctas
# set working directory
zctas = st_read("./data/Colorado_ZIP_Code_Tabulation_Areas_ZCTA.shp")
plot(st_geometry(zctas), axes = TRUE)
# convert longitude/latitude to UTM coordinates zone 13 using epsg
zctas = st_transform(zctas, 32613)
# make sure coordinate systems match up
plot(st_geometry(zctas), axes = TRUE)
points(st_coordinates(co), col = "orange", pch = 19)

# ensure coordinate reference systems are the same
st_crs(co) = st_crs(zctas)

# estimate omnidirectional semivariogram
# for aluminum (for simplicity, we won't
# examine anisotropy, but we should)
gco = gstat(id = "Al",
            formula = Al ~ ns(easting, knots = 5),
            data = co)

# crude fit
vco = variogram(gco, cutoff = 350000)
fit.vco = fit.variogram(vco,
                        vgm(psill = 0.5, model = "Sph",
                            range = 3e5, nugget = 1),
                        fit.method = 2)
plot(vco, fit.vco)

# create gstat object with estimated semivariogram
gco = gstat(id = "Al",
            formula = Al ~ ns(easting, knots = 5),
            data = co,
            model = fit.vco)

# point kriging prediction
centroids <- st_coordinates(st_centroid(zctas))
centroids <- st_as_sf(data.frame(easting = centroids[,1],
                                 northing = centroids[,2]),
                      coords = c("easting", "northing"),
                      remove = FALSE,
                      crs = st_crs(co))
# predict at centroids
p_point = predict(gco, newdata = centroids)

# block kriging
# add covariate columns
zctas$easting <- centroids$easting
zctas$northing <- centroids$northing
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