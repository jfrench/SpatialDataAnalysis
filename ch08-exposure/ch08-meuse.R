# install.packages(c("gstat", "geoR", "fields", "maptools", "gridExtra"))
library(sp)
library(gstat) # for most of the work

data(meuse, package = "sp") #available in sp package
# turn smoky dataframe into SpatialPointsDataFrame by adding coordinates
coordinates(meuse) <- c("x", "y")

# we'll do analysis on log(lead)
qqnorm(log(meuse$lead))

# fit variogram
v <- variogram(log(lead) ~ 1, meuse)
plot(v)
# model semivariogram with spherical semivariogram
f = fit.variogram(v, vgm(.5, "Sph", 1000, .1))
plot(v, f) # this fit looks goods

# this shows that we should really account for anisotropy, but we won't for simplicity
vd = variogram(log(lead) ~ 1, meuse, alpha = c(25, 70, 115, 160))
plot(vd)

# ordinary kriging on log(lead)
gmeuse = gstat(id = "lead", formula = log(lead) ~ 1,
               data = meuse, model = f)
data(meuse.grid) # prediction grid built into gstat
coordinates(meuse.grid) = ~ x + y
gridded(meuse.grid) = TRUE
# predict response
oklog = predict(gmeuse, newdata = meuse.grid)
spplot(oklog, "lead.pred")

# biased prediction, just take exp of predicted values
epred = exp(oklog$lead.pred)
oklog$lead.pred = epred
spplot(oklog, "lead.pred", col.regions = hcl.colors(64),
       cuts = 63)

# lognormal prediction
# lambda = 0 means use a log transform.  Related to box-cox transformation
olk = krigeTg(lead ~ 1, locations = meuse, newdata = meuse.grid,
              model = f, lambda = 0)
names(olk) # we want the trans-gaussian kriging predictions
spplot(olk, "var1TG.pred", col.regions = hcl.colors(64),
       cuts = 63)

# compare prediction maps on original scale (make sure colors on same scale)
cut = seq(0, 600, len = 63) # for consistent coloring of graphics
expplot = spplot(oklog, "lead.pred", colorkey = TRUE,
                 col.regions = hcl.colors(64), at = cut,
                 main = "biased predictions of lead")
olkplot = spplot(olk, "var1TG.pred", colorkey = TRUE,
                 col.regions = hcl.colors(64), at = cut,
                 main = "lognormal predictions of lead")
# see biasedness more clearly
head(cbind(oklog$lead.pred, olk$var1TG.pred))
head(cbind(oklog$lead.pred - olk$var1TG.pred))
library(gridExtra)
grid.arrange(expplot, olkplot, ncol = 2)


