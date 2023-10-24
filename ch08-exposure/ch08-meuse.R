# install.packages(c("gstat", "geoR", "fields", "maptools", "gridExtra"))
library(sp)
library(gstat) # for most of the work

data(meuse.all, package = "gstat") #available in sp package
# turn meuse.all dataframe into sf object
meuse <- sf::st_as_sf(meuse.all, coords = c("x", "y"))

# plot lead variable
plot(meuse["lead"], pal = hcl.colors, pch = 20)

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
meuse_grid <- st_as_sf(meuse.grid, coords = c("x", "y"))

# predict response
oklog = predict(gmeuse, newdata = meuse_grid)
plot(oklog["lead.pred"], pal = hcl.colors)

# biased prediction, just take exp of predicted values
epred = exp(oklog$lead.pred)
oklog$lead.pred = epred
plot(oklog["lead.pred"], pal = hcl.colors)

# lognormal prediction
# lambda = 0 means use a log transform.  Related to box-cox transformation
olk = krigeTg(lead ~ 1,
              locations = meuse,
              newdata = meuse_grid,
              model = f, lambda = 0)
names(olk) # we want the trans-gaussian kriging predictions
plot(olk["var1TG.pred"], pal = hcl.colors)

# join biased and unbiased predictions
Tgpred <- st_join(oklog, olk)

# plot results side by side
plot(Tgpred[c("lead.pred", "var1TG.pred")],
     pal = hcl.colors, key.pos = 4)

# see biasedness more clearly
head(cbind(oklog$lead.pred, olk$var1TG.pred))
head(cbind(oklog$lead.pred - olk$var1TG.pred))
