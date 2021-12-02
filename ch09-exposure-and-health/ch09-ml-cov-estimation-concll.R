library(mvtnorm) # for rmvnorm
library(gear) # for evgram

### covariance parameter estimation for non-stationary mean

set.seed(99)

N = 50

# observed data locations
coords = matrix(runif(N * 2, -1, 1), ncol = 2)

# distances between observed data locations
d = as.matrix(dist(coords))

# determine mean-related parameters
X = cbind(1, coords[,1])
beta = c(1, 2)

# construct isotropic exponential covariance function
# psill = partial sill
# a = range parameter
# c0 = nugget
# note: this function is incorrect if the coordinates are not unique
cov.exp = function(d, psill, a, c0) {
  psill * exp(-d/a) + c0 * diag(N)
}

S = cov.exp(d, psill = 1, a = 0.5, c0 = 0.1)

# generate observed data
y = c(rmvnorm(n = 1, mean = X %*% beta, sigma = S))

# negative concentrated multivariate normal log-likelihood
# thetastar = c(a, u), the covariance parameter
# vector.  needs to be a vector for nlminb function
conc.ll = function(thetastar) {
  a = thetastar[1]
  u = thetastar[2]

  # correlation matrix after dividing by psill
  V = cov.exp(d, psill = 1, a = a, c0 = 0) + u * diag(N)

  # estimate of beta (GLS)
  betatilde = solve(crossprod(X, solve(V, X)), crossprod(X, solve(V, y)))
  # residuals
  e = y - X %*% betatilde

  # estimate psill
  psilltilde = crossprod(e, solve(V, e))/N
  # -1 * concentrated log-likelihood of multivariate normal
  # -1 is because the nlminb function wants to minimize, not maximize
  1/2 * (N * log(2 * pi) +
         N * log(psilltilde) +
         determinant(V, log = TRUE)$mod +
         N)
}

# start is the starting parameter estimates
# objective is the function to minimize
# lower contrains the parameters to be no less than lower
# upper constrains the parameters to be no more than upper
nlminb(start = c(0.4, 0.01),
       objective = conc.ll,
       lower = c(0.01, 0),
       upper = c(2, 2))

# our "optimal" ML estimates of (a, u)
# $par
# [1] 0.2475395 0.0000000
# these are the maximum likelihood estimates of:
# a, c0/psill

# Thus, the MLE estimates are:
V = cov.exp(d, psill = 1, a = 0.2475395, c0 = 0) + 0 * diag(N)
(betatilde = solve(crossprod(X, solve(V, X)), crossprod(X, solve(V, y))))
# residuals
e = y - X %*% betatilde
(psilltilde = crossprod(e, solve(V, e))/N)
(c0 = 0/psilltilde)

# combined ml estimates
c(beta = betatilde, psill = psilltilde, a = 0.2475395, c0 = c0)

# same results using gear package
# formula: trend formula
# data: data frame with data
# coordnames: names of coordinate columns
# mod covariance model:
# r = range parameter (a)
# evar = c0 (error variance for filter model)
geolmod <- geolm(formula = y ~ x1, data = df, coordnames = c("x1", "x2"),
                 mod = cmod_std(model = "exponential", psill = 0.5,
                                r = 0.2, evar = 0.1))
# estimate parameters using ML (not REML)
geolmodtilde <- estimate(geolmod, reml = FALSE)
geolmodtilde$coeff
geolmodtilde$mod[c("psill", "r", "evar")]
