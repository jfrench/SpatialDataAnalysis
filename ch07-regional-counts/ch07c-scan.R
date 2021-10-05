library(spdep)
library(rgdal)
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

# observed cases, rounded down
cases = floor(ny8$Cases)
# population (same as nydf$population)
population = ny8$POP8
# expected number of cases
e = sum(cases)/sum(population) * population
# apply circular scan method
scan = scan.test(coords = coords,
                 cases = cases,
                 pop = population,
                 ex = e,
                 nsim = 999,
                 alpha  = 0.2)

# results from the test are available in
summary(scan)
# cluster information
clusters(scan)

# library(RColorBrewer) # useful for determining plotting colors
# # look at qualtitative color mapping that is colorblind friendly
# display.brewer.all(type = "qual", colorblindFriendly = TRUE)
# mycol = brewer.pal(3, "Dark2")
# # create vector of colors to show results
# # default is white (no clustering)
# nycol = rep("white", nrow(nydf))
# # the most likely cluster locations are lightblue
# nycol[scan$clusters[[1]]$locids] = mycol[1]
# # the secondary cluster locations are lightorange and lightpurple
# nycol[scan$clusters[[2]]$locids] = mycol[2]
# nycol[scan$clusters[[3]]$locids] = mycol[3]

# need to color 3 clusters
mycol = grDevices::hcl.colors(3)
# color.clusters(scan, col = mycol) colors the 3 clusters using the desired clusters
plot(sf::st_geometry(ny8), border="grey60", axes=TRUE,
     col = color.clusters(scan, col = mycol))
legend("topright", legend = c("Cluster A", "Cluster B", "Cluster C"),
       lwd = 10, col = mycol)

# a simpler plot
plot(scan)
