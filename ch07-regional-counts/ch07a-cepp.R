# install.packages("spdep", "rgdal", "RColorBrewer", "SpatialEpi)
library(spdep)
library(rgdal)
library(smerc)

# note that the Besag-Newell method in SpatialEpi is modified from the
# description in the book.
setwd("~/Dropbox/UCD_Files/Teaching/Math4027/Data/NY_data")
NY8 <- readOGR(".", "NY8_utm18")
NY_nb <- read.gal("NY_nb.gal", region.id=row.names(NY8))

# plot regions
plot(NY8, border="grey60", axes=TRUE)
# plot neighbors
plot(NY_nb, coordinates(NY8), pch = 19, cex = 0.6, add = TRUE)

# read in data
setwd("~/Dropbox/UCD_Files/Teaching/Math4027/Data/")
nydf = read.table("NYtract.dat")
names(nydf) = c("x", "y", "Population", "Observed")
coords = nydf[,c("x", "y")]

# note that is you use the x, y coordinates of NY8 (coordinates(NY8))
# you may get different results since the coordinates
# have been projected to the UTM coordinate system
# We have sought to match the book's results

# Use CEPP method for different values of n*
# coords is the x,y coordinates
# pop is population size
# cases is the number of cases
# nstar is the number of at-risk person to include in the circle radius
# alpha is the significance level for which to identify a cluster
cepp1000 = cepp.test(coords = coords,
              cases = nydf$Observed,
              pop = nydf$Population,
              nstar = 1000,
              alpha = 0.10)
cepp5000 = cepp.test(coords = coords,
                     cases = nydf$Observed,
                     pop = nydf$Population,
                     nstar = 5000,
                     alpha = 0.10)

cepp10000 = cepp.test(coords = coords,
                     cases = nydf$Observed,
                     pop = nydf$Population,
                     nstar = 10000,
                     alpha = 0.10)

cepp40000 = cepp.test(coords = coords,
                     cases = nydf$Observed,
                     pop = nydf$Population,
                     nstar = 40000,
                     alpha = 0.10)

# Note:  we generally get a lot of clusters
# we only look at the most likely one
# look at x$clusters to see them all

# plot likeliest clusters
plot(NY8, border = "grey60", axes = TRUE, col = color.clusters(cepp1000))
legend("topright", legend = c("n* = 1000"))
plot(cepp1000)

plot(NY8, border = "grey60", axes = TRUE, col = color.clusters(cepp5000))
legend("topright", legend = c("n* = 5000"))
plot(cepp5000)

plot(NY8, border="grey60", axes = TRUE, col = color.clusters(cepp10000))
legend("topright", legend = c("n* = 10000"))
plot(cepp10000)

plot(NY8, border="grey60", axes = TRUE, col = color.clusters(cepp40000))
legend("topright", legend = c("n* = 40000"))
plot(cepp40000)