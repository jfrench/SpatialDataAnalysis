states = map("state")
SpatialPolygonsDataFrame(states)

x = readShapePoly(file.choose())
nbx = poly2nb(x)
plot(nbx, coords = coordinates(x), add = TRUE)
