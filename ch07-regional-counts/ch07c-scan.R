# install.packages("spdep", "rgdal", "RColorBrewer", "smerc")
library(spdep)
library(rgdal)
library(smerc)

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

# note that is you use the x, y coordinates of NY8 (coordinates(NY8))
# you will get different results since the coordinates
# have been projected to the UTM coordinate system
# We have sought to match the book's results
geolonlat = nydf[,c("x", "y")]
# geo = latlong2grid(geolonlat)
# get cases, population, and expected numbers
cases = floor(NY8$Cases)
population = NY8$POP8
e = sum(cases)/sum(population)*population
scan = scan.test(coords = geolonlat,
                 cases = cases, 
                 pop = population, 
                 ex = e,
                 nsim = 999,
                 alpha  = 0.2)

# results from the test are avaialble in 
scan$clusters

library(RColorBrewer) # useful for determining plotting colors
# look at qualtitative color mapping that is colorblind friendly
display.brewer.all(type = "qual", colorblindFriendly = TRUE)
mycol = brewer.pal(3, "Dark2")
# create vector of colors to show results
# default is white (no clustering)
nycol = rep("white", nrow(nydf))
# the most likely cluster locations are lightblue
nycol[scan$clusters[[1]]$locids] = mycol[1]
# the secondary cluster locations are lightorange and lightpurple
nycol[scan$clusters[[2]]$locids] = mycol[2]
nycol[scan$clusters[[3]]$locids] = mycol[3]

# plot most likely clusters
plot(NY8, border="grey60", axes=TRUE, col = nycol)
legend("topright", legend = c("Cluster A", "Cluster B", "Cluster C"), 
       lwd = 10, col = mycol)