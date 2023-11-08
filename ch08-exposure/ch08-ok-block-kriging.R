# install.packages(c("gstat", "geoR", "fields", "maptools", "gridExtra"))
library(gstat) # for most of the work
library(sp) # for Spatial... obects
library(rgdal) # to read shapefile and spTransform
library(colorspace) # for terrain_hcl

# load colorado data
data(co, package = "gear")
coordinates(co) = c("easting", "northing")
spplot(co, "Al")

# read colorado zctas
# set working directory
zctas = readOGR("./data/co_zcta", "Colorado_ZIP_Code_Tabulation_Areas_ZCTA")
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

# point kriging prediction
coordnames(zctas) = c("easting", "northing")
# have to make coordinate names match
coords = coordinates(zctas)
colnames(coords) = c("easting", "northing")
centroids = SpatialPoints(coords,
                          proj4string = CRS(proj4string(zctas)))
p_point = predict(gco, newdata = centroids)

# block kriging
p_block = predict(gco, newdata = zctas)

# plot of block kriging
spplot(p_block, "Al.pred", cuts = 11,
       col.regions = colorspace::terrain_hcl(12))

# compare to point prediction
spplot(p_point, "Al.pred", cuts = 11,
       col.regions = colorspace::terrain_hcl(12),
       key.space = "right")
