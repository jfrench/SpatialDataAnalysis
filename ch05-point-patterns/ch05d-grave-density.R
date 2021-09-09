library(spatstat) # for inference on spatial point processes
library(smacpod)
data(grave, package = "smacpod") # import data in ppp format

# determine affected and unaffected sides
af <- which(grave$marks == "affected")
un <- which(grave$marks == "unaffected")

# plot of event locations
plot(grave)

# recommended bandwidths for affected and unaffected graves
# in u- and v-directions
dim <- 2
# Scott's bandwidths in u- and v-directions
# for each group
buaf <- sd(grave$x[af])*length(af)^(-1/(dim+4))
bvaf <- sd(grave$y[af])*length(af)^(-1/(dim+4))
buun <- sd(grave$x[un])*length(un)^(-1/(dim+4))
bvun <- sd(grave$y[un])*length(un)^(-1/(dim+4))

# density estimated for affected and unaffected groups
# using associated bandwidths from Scott's rule
iaf <-  spdensity(grave[af,], sigma = c(buaf, bvaf))
iun <-  spdensity(grave[un,], sigma = c(buun, bvun))

# plot perspective and contour plots of estimated density for affected sites
par(mfrow = c(1, 2))
persp(iaf, theta = 45, phi = 35, xlab = "u", ylab = "v", zlab = "density", main = "Estimated intensity function")
contour(iaf, xlab = "u", ylab = "v", main = "Affected grave locations")
points(grave, pch = ".")

# plot perspective and contour plots of estimated density for unaffected sites
par(mfrow = c(1, 2))
persp(iun, theta = 45, phi = 35, xlab = "u", ylab = "v", zlab = "density", main = "Estimated intensity function")
contour(iun, xlab = "u", ylab = "v", main = "Unaffected grave locations")
points(grave, pch = ".")

# plot perspective and contour plots of estimated density for unaffected sites
par(mfrow = c(1, 2))
contour(iaf, xlab = "u", ylab = "v", main = "Affected")
points(grave, pch = ".")
contour(iun, xlab = "u", ylab = "v", main = "Unaffected")
points(grave, pch = ".")

par(mfrow = c(1, 1))
