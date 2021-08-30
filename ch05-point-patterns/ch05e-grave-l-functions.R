library(spatstat) # for inference on spatial point processes
# source("http://math.ucdenver.edu/~jfrench/data/ss.R") # for custom helper functions

setwd("~/OneDrive - The University of Colorado Denver/teaching/math6384/sda_code/data")
grave = read.csv("grave.csv", header = TRUE)
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
# the tolerance envelopes
# level is the tolerance level of the envelopes
# ... is additional arguments passed to Lest
lplot <- function(x, nsim = 500, level = 0.95,
                  correction = "Ripley", ...) {
  lobs <- spatstat.core::Lest(x, correction = correction, ...)

  # lambda <- summary(x)$intensity
  win <- x$window
  lsim <- pbapply::pblapply(1:nsim, FUN = function(i) {
    # xsim <- spatstat::rpoispp(lambda, win = win)
    # generate n event locations over study area under CSR
    xsim <- spatstat.core::runifpoint(n = x$n, win = win)
    # estimate L for simulated point pattern
    spatstat.core::Lest(xsim, correction = correction, ...)
  })

  r <- lobs$r # get distances
  obs <- lobs$iso # get estimated l for observed
  # get estimated l for each simulated data set
  sim <- sapply(lsim, getElement, "iso")
  # apply the min function to each row  (MARGIN = 1) of sim
  # gets pointwise minimum for simulated data
  # at each distance.  do same for max, quantiles, median
  lo <- apply(sim, MARGIN = 1, FUN = min, na.rm = TRUE)
  hi <- apply(sim, MARGIN = 1, FUN = max, na.rm = TRUE)
  alpha <- 1 - level
  qlo <- apply(sim, MARGIN = 1, FUN = quantile,
               prob = alpha/2, na.rm = TRUE)
  qhi <- apply(sim, MARGIN = 1, FUN = quantile,
               prob = 1 - alpha/2, na.rm = TRUE)
  med <- apply(sim, MARGIN = 1, FUN = median, na.rm = TRUE)
  # construct empty plot of the right size
  plot(range(r), c(min(c(lo, obs) - r, na.rm = TRUE),
                   max(c(hi, obs) - r, na.rm = TRUE)),
       type = "n",
       xlab = "distance", ylab = "L(distance) - distance")
  # plot different statistics with different styles/thickness
  lines(r, obs - r, lwd = 2)
  lines(r, lo - r, lty = 2)
  lines(r, hi - r, lty = 2)
  lines(r, qlo - r, lty = 1)
  lines(r, qhi - r, lty = 1)
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

# Are the affected graves more clustered than expected
# if the grave sites were randomly labeled?

# estimate L - r for observed
# use Ripley's correction
laf = spatstat.core::Lest(pppol[,af], r = r, correction = "Ripley")$iso - r # Lhat for affected

# get n and naf for simulations
n <- pppol$n
naf <- length(af)

# 499 times, randomly label affected locations and compute
# L - r for affected sites
# pbsapply returns this in a nice format and includes a timer
# each column is a simulation
lrl <- pbapply::pbsapply(1:499, FUN = function(i) {
  xsim <-
  spatstat.core::Lest(pppol[,sample.int(n, size = naf)],
       r = r, correction = "Ripley")$iso - r
})

# apply the min function to each row  (MARGIN = 1) of lrl
lo <- apply(lrl, MARGIN = 1, FUN = min, na.rm = TRUE)
hi <- apply(lrl, MARGIN = 1, FUN = max, na.rm = TRUE)
qlo <- apply(lrl, MARGIN = 1, FUN = quantile,
             prob = 0.025, na.rm = TRUE)
qhi <- apply(lrl, MARGIN = 1, FUN = quantile,
             prob = 0.975, na.rm = TRUE)
med <- apply(lrl, MARGIN = 1, FUN = median, na.rm = TRUE)
# construct empty plot of the right size
plot(c(0, 5000), c(-600, 800), type = "n",
     xlab = "distance", ylab = "L(distance) - distance")
lines(r, laf, lwd = 2)
lines(r, lo, lty = 2)
lines(r, hi, lty = 2)
lines(r, qlo, lty = 1)
lines(r, qhi, lty = 1)
lines(r, med, lty = 4)
legend("topleft",
       legend = c("min/max",
                  "2.5th, 97.5th percentiles",
                  "median"),
       lty = c(2, 1, 4))

# test hypothesis that affected event locations are
# clustered at any scale for 0 <= h <= 2000
# in rectangular domain
r <- seq(0, 2000, len = 201)
Tobs <- max(spatat.core::Lest(pppafpol, r = r, correction = "Ripley")$iso - r)
# relabel affected events,
# then compute max(Lhat(h) - h)
# for relabeled data
Tsim <- pbapply::pbsapply(1:499, FUN = function(i) {
  max(spatstat.core::Lest(pppol[sample.int(n, size = naf), 2:3], r = r,
           correction = "Ripley")$iso - r)
})

# proportion of simulated test statistics
# as extreme as one observed
# the observed pattern is relative consistent with a random labeling of affected
mean(c(Tsim, Tobs) >= Tobs)
