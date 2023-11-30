library(mvtnorm) # for rmvnorm, dmvnorm
library(gear) # for evgram

### covariance parameter estimation for non-stationary mean

set.seed(99)

# observed data locations
coords = matrix(runif(100, -1, 1), ncol = 2)

# distances between observed data locations
d = as.matrix(dist(coords))

# determine mean-related parameters
X = cbind(1, coords[,1])
beta = c(1, 2)

# construct isotropic exponential covariance function
# note: this function is incorrect if the coordinates are not unique
cov.exp = function(d, psill, range, nugget) {
  psill * exp(-d/range) + nugget * diag(nrow(d))
}

S = cov.exp(d, psill = 1, range = 0.5, nugget = 0.1)

# generate observed data
y = c(rmvnorm(n = 1, mean = X %*% beta, sigma = S))

# Negative concentrated multivariate normal log-likelihood (only beta profiled)
# theta = c(psill, a, c0), the covariance parameter vector.
# Needs to be a vector for nlminb function
conc.ll = function(theta) {
  psill = theta[1]
  a = theta[2]
  c0 = theta[3]

  # covariance matrix based on given parameters
  Shat = cov.exp(d, psill, a, c0)
  # estimate of beta (GLS)
  betahat = solve(crossprod(X, solve(Shat, X)),
                  crossprod(X, solve(Shat, y)))
  # log-likelihood of multivariate normal using mvtnorm
  # - is because the nlminb function wants to minimize,
  # not maximize
  - dmvnorm(y, mean = X %*% betahat, Shat, log = TRUE)
}

# construct empirical semivariogram using gear package to
# find reasonable starting values for covariance parameters
df = data.frame(y, x1 = coords[,1], x2 = coords[,2])
# detrend data before estimating semivariogram
ev = evgram(y ~ x1, data = df, coordnames = ~ x1 + x2)
plot(ev, ylim = c(0, 0.7))

# start is the starting parameter estimates
# objective is the function to minimize
# lower constrains the parameters to be no less than lower
# upper constrains the parameters to be no more than upper
nlminb(start = c(0.4, 0.4, 0.1),
       objective = conc.ll,
       lower = c(0.01, 0.01, 0),
       upper = c(2, 2, 2))

# our "optimal" ML estimates of the covariance are
# $par
# [1] 0.45118 0.24754 0.00000
# these are the maximum likelihood estimates of:
# partial sill, range parameter, and nugget
# covariance matrix based on given parameters
Shat = cov.exp(d, 0.45118, 0.24754, 0)
# estimate of beta (GLS)
(betatilde = solve(crossprod(X, solve(Shat, X)), crossprod(X, solve(Shat, y))))



