library(spatstat) # for inference on spatial point processes
# source("http://math.ucdenver.edu/~jfrench/data/ss.R") # for custom helper functions

grave = read.csv("./data/grave.csv", header = TRUE)
grave$deformity = factor(grave$deformity)
levels(grave$deformity) = c("unaffected", "affected")
# which event locations are affected/unaffected
af = which(grave$deformity == "affected")
un = which(grave$deformity == "unaffected")

#rec = locator()
# points must be anticlockwise!
# create boundary of rectangular window
delta <- 100
umin = min(grave[,2]) - delta; umax = max(grave[,2]) + delta
vmin = min(grave[,3]) - delta; vmax = max(grave[,3]) + delta
rec = data.frame(x = c(umax, umin, umin, umax), y = c(vmax, vmax, vmin, vmin))
# create boundary of polygonal window
#poly = locator()
poly = data.frame(x = c(9464.853,  8735.634,  8006.414,  8006.414,  6629.000,  4927.488,  4157.756,  4279.293,  4846.463,  5980.805,  6466.951,  7884.878,  7844.366,  8046.927,  8046.927,  9910.487, 10923.292, 10315.609, 10882.780, 10356.122, 10599.195,  9951.000, 10275.097),
                   y = c(10749.814,  9939.571,  9858.546,  9250.863,  9088.814,  7751.912,  7711.400,  6090.912,  5523.742,  5564.254,  6212.449,  5321.181,  4591.961,  4429.912,  3943.766,  2768.912, 3093.010,  3660.181,  4227.351,  5604.766,  8035.497,  8764.717,  9372.400))

wrec = owin(poly = rec) # create rectangular bounding window
wpoly = owin(poly = poly) # create polygon bounding window

# create point pattern object with rectangular and polygonal boundaries
pprec = ppp(x = grave[,2], y = grave[,3], window = wrec, marks = factor(grave[,1]))
pppol  = ppp(x = grave[,2], y = grave[,3], window = wpoly, marks = factor(grave[,1]))
pppunrec = as.ppp(grave[un,2:3], wrec)
pppafrec = as.ppp(grave[af,2:3], wrec)
pppunpol = as.ppp(grave[un,2:3], wpoly)
pppafpol = as.ppp(grave[af,2:3], wpoly)

# distances to evaluate K function
r = seq(0, 5200, len = 512)
set.seed(1)

# x is a ppp
# nsim is the number of simulation from which to construct
# the envelopes
# level is the quantile level of the envelopes
# ... is additional arguments passed to Lest
lplot <- function(x, nsim = 499, level = 0.95,
                  correction = "Ripley", ...) {
  e_outer = envelope(x, fun = spatstat.explore::Lest,
                     nsim = nsim, nrank = 1,
                     savepatterns = TRUE)
  e_tol = envelope(e_outer, fun = spatstat.explore::Lest,
                   nsim = nsim,
                   nrank = floor((1 - level) * (nsim + 1)/2))
  plot(e_outer, fmla = . - r ~ r, legend = FALSE)
  plot(e_tol, fmla = . - r ~ r, shadecol = "lightgrey", add = TRUE)
}

## Are event locations clustered in the rectangle?
## Fix the number of events in the study area
# L plot for rectangle
lplot(pprec, r = r)
title("L plot for all grave sites, rectangle")

# L plot for affected/rectangle
lplot(pppafrec, r = r)
title("L plot for affected grave sites, rectangle")

# L plot for unaffected/rectangle
lplot(pppunrec, r = r)
title("L plot for nonaffected grave sites, rectangle")

# Are sites clustered within a different polygon?
# L plots for polygon
lplot(pppol, r = r)
title("L plot for all grave sites, polygon")

# L plots for affected/polygon
lplot(pppafpol, r = r)
title("L plot for affected grave sites, polygon")

# L plots for unaffected/polygon
lplot(pppunpol, r = r)
title("L plot for nonaffected grave sites, polygon")

# test hypothesis that affected event locations are
# clustered at any scale for 0 <= h <= 2000
# in polygon domain
r <- seq(0, 2000, len = 201)
Tobs <- max(spatstat.explore::Lest(pppafpol, r = r, correction = "Ripley")$iso - r)
# relabel affected events,
# then compute max(Lhat(h) - h)
# for relabeled data
Tsim <- pbapply::pbsapply(1:499, FUN = function(i) {
  max(spatstat.explore::Lest(pppol[sample.int(n, size = naf), 2:3], r = r,
           correction = "Ripley")$iso - r)
})

# proportion of simulated test statistics
# as extreme as one observed
# the observed pattern is relative consistent with a random
# labeling of affected
mean(c(Tsim, Tobs) >= Tobs)
