library(HistData) # load package with John Snow cholera data

data(Snow.polygons) # polygons surrounding each water pump

# determine range of x and y values to know how large the plot should be
extracted_x <- sapply(Snow.polygons, getElement, "x")
range_x <- range(unlist(extracted_x))
extracted_y <- sapply(Snow.polygons, getElement, "y")
range_y <- range(unlist(extracted_y))

# create an empty plot
plot(range_x, range_y, type = "n", xlim = range_x, ylim = range_y)
lapply(Snow.polygons, function(df) {
  polygon(df$x, df$y)
})

# load data set with locations of cholera death
data(Snow.deaths)
# plot cholera death locations
points(Snow.deaths$x, Snow.deaths$y, pch = 19, col = "lightgrey")

# load data set with water pump locations
data(Snow.pumps)
# plot water pump locations
points(Snow.pumps$x, Snow.pumps$y, pch = 4, col = "blue", cex = 3)
