# install.packages("spdep", "sf", "smerc", "RColorBrewer")
library(spdep)
library(sf)
library(smerc)

# read shapefile for new york counties
ny8 <- sf::st_read("./data/NY_data/NY8_utm18.shp")
# read neighbor information
ny_nb <- spdep::read.gal("./data/NY_data/NY_nb.gal", override.id = TRUE)

# plot region boundaries from ny8
plot(st_geometry(ny8), border="grey60")
# plot neighbors
plot(ny_nb, coords = st_centroid(st_geometry(ny8)),
     add=TRUE, col="blue", pch = 19, cex = 0.6)

# read in data
nydf = read.table("./data/NYtract.dat")
names(nydf) = c("x", "y", "Population", "Observed")
coords = nydf[,c("x", "y")]

# note that is you use the x, y coordinates of ny8 (sf::st_centroids(ny8))
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
# pretty plot
plot(sf::st_geometry(ny8), border = "grey60", axes = TRUE,
     col = color.clusters(cepp1000))
legend("topright", legend = c("n* = 1000"))
# basic plot
plot(cepp1000)
# basic info
cepp1000
# cluster info
summary(cepp1000)
clusters(cepp1000)


plot(sf::st_geometry(ny8), border = "grey60", axes = TRUE,
     col = color.clusters(cepp5000))
legend("topright", legend = c("n* = 5000"))
plot(cepp5000)
# basic info
cepp5000
# cluster info
summary(cepp5000)
clusters(cepp5000)

plot(sf::st_geometry(ny8), border="grey60", axes = TRUE,
     col = color.clusters(cepp10000))
legend("topright", legend = c("n* = 10000"))
plot(cepp10000)
# basic info
cepp10000
# cluster info
summary(cepp10000)
clusters(cepp10000)

plot(sf::st_geometry(ny8), border="grey60", axes = TRUE,
     col = color.clusters(cepp40000))
legend("topright", legend = c("n* = 40000"))
# basic info
cepp40000
# cluster info
summary(cepp40000)
clusters(cepp40000)