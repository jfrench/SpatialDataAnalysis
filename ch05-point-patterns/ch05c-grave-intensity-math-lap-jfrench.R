library(spatstat) # for inference on spatial point processes

# import and clean data
setwd("./data/")
grave <- read.csv("grave.csv", header = TRUE)
grave$deformity <-  factor(grave$deformity)
levels(grave$deformity) <-  c("unaffected", "affected")

# determine affected and unaffected sides
af <- which(grave$deformity == "affected")
un <-  which(grave$deformity == "unaffected")

# plot of event locations
plot(v ~ u, data = grave)

# rec = locator()
# points must be anticlockwise!
# create boundary of rectangular window
umin <-  min(grave[,2]) - 500; umax = max(grave[,2]) + 500
vmin <-  min(grave[,3]) - 500; vmax = max(grave[,3]) + 500
# create rectangular bounding window
rec <-  data.frame(x = c(umax, umin, umin, umax),
                 y = c(vmax, vmax, vmin, vmin))
wrec <-  owin(poly = rec) # create rectangular bounding window
# create point pattern object with rectangular and polygonal boundaries
pprec <-  ppp(x = grave[,2], y = grave[,3],
            window = wrec, marks = factor(grave[,1]))

# plot affected and unaffected individuals in study area
plot(pprec, main = "grave locations")

# recommended bandwidths for affected and unaffected graves
# in u- and v-directions
dim <-  2
# Scott's bandwidths in u- and v-directions
# for each group
buaf <-  sd(grave[af,2])*length(af)^(-1/(dim+4))
bvaf <-  sd(grave[af,3])*length(af)^(-1/(dim+4))
buun <-  sd(grave[un,2])*length(un)^(-1/(dim+4))
bvun <-  sd(grave[un,3])*length(un)^(-1/(dim+4))

# intensity estimated for affected and unaffected groups
# using associated bandwidths from Scott's rule
iaf <-  density(pprec[af,], sigma = c(buaf, bvaf))
iun <-  density(pprec[un,], sigma = c(buun, bvun))

# plot perspective and contour plots of estimated intensity for affected sites
par(mfrow = c(1, 2))
persp(iaf, theta = 45, phi = 35, xlab = "u", ylab = "v", zlab = "density", main = "Estimated intensity function")
contour(iaf, xlab = "u", ylab = "v", main = "Affected grave locations")
points(grave[af,2:3], pch = ".")

# plot perspective and contour plots of estimated intensity for unaffected sites
par(mfrow = c(1, 2))
persp(iun, theta = 45, phi = 35, xlab = "u", ylab = "v", zlab = "density", main = "Estimated intensity function")
contour(iun, xlab = "u", ylab = "v", main = "Unaffected grave locations")
points(grave[un, 2:3], pch = ".")

# plot perspective and contour plots of estimated intensity for affecteed/unaffected sites jointly
par(mfrow = c(1, 2))
contour(iaf, xlab = "u", ylab = "v", main = "Affected")
points(grave[af,2:3], pch = ".")
contour(iun, xlab = "u", ylab = "v", main = "Unaffected")
points(grave[un, 2:3], pch = ".")

par(mfrow = c(1, 1)) # set back to 1x1 plot

# manually reproduce a plot
kern <-


