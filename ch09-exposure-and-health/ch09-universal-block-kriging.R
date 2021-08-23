# install.packages(c("gstat", "geoR", "fields", "maptools", "gridExtra"))
library(gstat) # for most of the work
library(sp) # for Spatial... obects
library(rgdal) # to read shapefile and spTransform
library(splines) # for ns
library(colorspace) # for terrain_hcl

# load colorado data
data(co, package = "gear")
coordinates(co) = c("easting", "northing")
spplot(co, "Al")

# read colorado zctas
# set working directory
setwd("~/OneDrive - The University of Colorado Denver/Teaching/Math6384/Data/co_zcta")
zctas = readOGR(".", "Colorado_ZIP_Code_Tabulation_Areas_ZCTA")
plot(zctas, axes = TRUE)
# convert longitude/latitude to UTM coordinates zone 13
zctas = spTransform(zctas, "+proj=utm +zone=13 +datum=NAD83")
plot(zctas, axes = TRUE)

# ensure coordinate reference systems are the same
proj4string(co) = CRS(proj4string(zctas))

# estimate omnidirectional semivariogram
# for aluminum (for simplicity, we won't
# examine anisotropy, but we should)
gco = gstat(id = "Al",
            formula = Al ~ ns(easting, knots = 5),
            data = co)

# crude fit
vco = variogram(gco, cutoff = 350000)
fit.vco = fit.variogram(vco, vgm(psill = 0.5, model = "Sph",
                                 range = 3e5, nugget = 1),
                        fit.method = 2)
plot(vco, fit.vco)

# create gstat object with estimated semivariogram
gco = gstat(id = "Al",
            formula = Al ~ ns(easting, knots = 5),
            data = co,
            model = fit.vco)

# point kriging prediction
coordnames(zctas) = c("easting", "northing")
# have to make coordinate names match
coords = coordinates(zctas)
colnames(coords) = c("easting", "northing")
centroids = SpatialPoints(coords,
                          proj4string = CRS(proj4string(zctas)))
p_point = predict(gco, newdata = centroids)


# block kriging prediction
# create vector for grid for each polygon
poly_grids = vector("list", length(zctas@polygons))
for (i in seq_along(zctas@polygons)) {
  # predict on grid within each region
  border = zctas@polygons[[1]]@Polygons[[1]]
  g = spsample(zctas[i,], type = "regular", n = 50)
  coordnames(g) = c("easting", "northing")
  poly_grids[[i]] = g
}

# number of points in grid for each zcta
ng = sapply(poly_grids, function(x) nrow(x@coords))
# collapse grid points
newdata = sapply(poly_grids, function(x) x@coords)
# rbind grid points for each component of the list
newdata = do.call(rbind, poly_grids)

# takes awhile since there are a lot of points
p_block = predict(gco, newdata = newdata)
# average predictions for each polygon
avg_Al = tapply(X = p_block$Al.pred,
                INDEX = rep(seq_along(ng), times = ng),
                FUN = mean)
# add average value to zctas SpatialPolygonsDataFrame
zctas$avg_Al = avg_Al
# plot results for block kriging
spplot(zctas, "avg_Al", cuts = 11,
       col.regions = colorspace::terrain_hcl(12))
# compare to point prediction
spplot(p_point, "Al.pred", cuts = 11,
       col.regions = colorspace::terrain_hcl(12),
       key.space = "right")
