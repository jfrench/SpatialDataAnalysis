set.seed(1)

# load packages
library(smacpod)
data(grave)

# selection to graves that are affected
af = which(grave$marks == "affected")
# estimate densities of cases and controls, respectively,
# using bandwidth of 700
f700 = spdensity(grave[af, ], sigma = 700)
g700 = spdensity(grave[-af,], sigma = 700)

# contour plot of f700 with title
contour(f700, nlevels = 15, main = "")
title("Affected, bandwidth = 700")
# contour plot of g700, with title
contour(g700, nlevels = 15, main = "")
title("Nonaffected, bandwidth = 700")
# log ratio of spatial densities
r700 = logrr(grave, case = "affected", sigma = 700)
# contour plot of r700 (lty and lwd determined
# by experimentation)
contour(r700, lty = c(1, 1, 1, 1, 2, 1),
        lwd = c(1, 1, 1, 1, 1, 2), main = "")
title("Gaussian kernel, Bandwidth = 700")

# calculate log ratio of spatial densities for bandwidth = 350
f350 = spdensity(grave[af, ], sigma = 350)
g350 = spdensity(grave[-af,], sigma = 350)
r350 = logrr(grave, case = "affected", sigma = 350)
# contour plot of log ratio of spatial densities
contour(r350, lty = c(2, 1, 1, 1, 1), main = "")
# construct 95% tolerance envelopes for log relative risk
# when bandwidth = 350
if (!file.exists("renv350.rda")) {
  renv350 = logrr(grave, sigma = 350, case = "affected",
                  nsim = 999, level = 0.95)
  save(renv350, file = "renv350.rda")
}
load("renv350.rda")
# image plot showing regions where r350 is outside
# tolerance envelopes
plot(renv350)

# a better color scale (in my opinion)
# making it easier to distguish the clusters of cases relative
# to controls (red) and vice versa (blue)
grad = gradient.color.scale(min(renv350$v, na.rm = TRUE),
                            max(renv350$v, na.rm = TRUE))
plot(renv350, col = grad$col, breaks = grad$breaks)

## global test that spatial densities of cases/controls
# are the same
logrr.test(renv350)

# construct 95% tolerance envelopes for bandwidth of 700
if (!file.exists("renv700.rda")) {
  renv700 = logrr(grave, case = "affected", nsim = 999,
                  sigma = 700, nyx = c(60, 60),
                  level = 0.95)
  save(renv700, file = "renv700.rda")
}
load("renv700.rda")
# image plot showing regions where r700 is outside tolerance
# envelopes
plot(renv700)

### difference in K functions
# estimate K function for affected and unaffected,
# take their difference
kd = kdest(grave, case = "affected")
# plot estimated and theoretical KD (under RLH)
plot(kd, cbind(iso, theo) ~ r, legend = FALSE, main = "")

# construct envelopes using 499 randomly labeled data sets
# r is chosen to match the book example. In general,
# you shouldn't specify r unless you know what spatial scales
# you want to consider.
nsim = 499
if (!file.exists("kdenv.rda")) {
  kdenv = kdest(grave, case = "affected", nsim = 499,
                r = seq(0, 2000, len = 201),
                level = 0.95)
  save(kdenv, file = "kdenv.rda")
}
load("kdenv.rda")
# print some information about kdenv
print(kdenv)
# determine distances KD(r) outside 95% tolerance envelopes
summary(kdenv)
# plot results
plot(kdenv, ylab = "difference in K functions",
     xlim = c(0, 2000))
# legend("topleft",
#        legend = c("obs", "avg", "max/min env", "95% env"),
#        lty = c(1, 2, 1, 2),
#        col = c("black", "red", "darkgrey", "lightgrey"),
#        lwd = c(1, 1, 10, 10))

# KD+ global test
# H0: KD(r) = 0 for all r considered
# Ha: KD(r) > 0 for at least one r considered
kdplus.test(kdenv)
