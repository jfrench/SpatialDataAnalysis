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

# Use besag-newell method for different values
# of cstar.  cstar = 6, 12, 17, 23
# coords is the x,y coordinates
# population size
# cases is the number of cases
# cstar is the number of cases to include in the circle radius
# alpha is the significance level for which to identify a cluster
bn6 = bn.test(coords = coords,
              cases = nydf$Observed,
              pop = nydf$Population,
              cstar = 6,
              alpha = 0.01)
bn12 = bn.test(coords = coords,
               cases = nydf$Observed,
               pop = nydf$Population,
               cstar = 12,
               alpha = 0.01)
bn17 = bn.test(coords = coords,
               cases = nydf$Observed,
               pop = nydf$Population,
               cstar = 17,
               alpha = 0.01)
bn23 = bn.test(coords = coords,
               cases = nydf$Observed,
               pop = nydf$Population,
               cstar = 23,
               alpha = 0.01)

# Note:  we generally get a lot of clusters
# we only look at the most likely one
# look at x$clusters to see them all

# notice that the most likely clusters for different k
# largely overlap
bn6$clusters[[1]]$locids
bn12$clusters[[1]]$locids
bn17$clusters[[1]]$locids
bn23$clusters[[1]]$locids

library(RColorBrewer) # useful for determining plotting colors
# look at qualitative color mapping that is colorblind friendly
display.brewer.all(type = "qual", colorblindFriendly = TRUE)
mycol = brewer.pal(4, "Dark2")
# create vector of colors to show results
# default is white (no clustering)
nycol = rep("black", nrow(nydf))
# the most likely cluster locations are lightorange for cstar = 12
nycol[bn12$clusters[[1]]$locids] = mycol[2]
# the most likely cluster locations are lightgreen for cstar = 17
nycol[bn6$clusters[[1]]$locid] = mycol[1]
# the most likely cluster locations are magenta for cstar = 6, 17
nycol[bn6$clusters[[1]]$locid] = mycol[4]
# the most likely cluster locations are lightpurple for cstar = 23
nycol[bn23$clusters[[1]]$location] = mycol[3]

# plot most likely clusters
plot(NY8, border="grey60", axes = TRUE, col = nycol)
legend("topright",
       legend = c("Cluster k = 6, 17", "Cluster k = 12",
                  "Cluster k = 17", "Cluster k = 23"),
       lwd = 10, col = mycol)

# look more closely at most likely cluster information
bn6$clusters[[1]]
bn12$clusters[[1]]
bn17$clusters[[1]]
bn23$clusters[[1]]

# summary of results
# (the book seems to have incorrect p-values based on using
# c* instead of c* - 1 when calculating p-values)
# c*    observed    n_i,c*    local p-value
# 6     8.18        2921      0.007
# 12    12.21       8876      0.005
# 17    17.69       11268     0.0003
# 23    25.46       19615     0.0001
