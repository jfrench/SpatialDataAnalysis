library(mvtnorm)
library(gear)

### covariance parameter estimation for non-stationary mean

set.seed(99)

# observed data locations
coords = matrix(runif(100, -1, 1), ncol = 2)

# distances between observed data locations
d = as.matrix(dist(coords))

# determine mean-related parameters
X = cbind(1, coords[,1])
beta = c(1, 2)

# construct exponential covariance function
# note: this will not work for filtered kriging settings or
# when the data locations aren't unique
cov.exp = function(d, psill, range, nugget) {
  psill * exp(-d/range) + nugget * diag(nrow(d))
}

S = cov.exp(d, psill = 1, range = 0.5, nugget = 0.1)

# generate observed data
y = c(rmvnorm(n = 1, mean = X %*% beta, sigma = S))

# negative concentrated multivariate normal log-likelihood
# theta = c(psill, range, nugget), the covariance parameter
# vector.  needs to be a vector for nlminb function
conc.ll = function(theta) {
  psill = theta[1]
  range = theta[2]
  nugget = theta[3]

  # covariance matrix based on given parameters
  Shat = cov.exp(d, psill, range, nugget)
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
ev = vgram(y ~ x1, data = df, coords = ~ x1 + x2)
plot(ev, ylim = c(0, 0.7))

# start is the starting parameter estimates
# objective is the function to minimize
# lower contrains the parameters to be no less than lower
# upper constrains the parameters to be no more than upper
nlminb(start = c(0.4, 0.4, 0.1),
       objective = conc.ll,
       lower = c(0.01, 0.01, 0),
       upper = c(2, 2, 2))

# our "optimal" ML estimates of the covariance are
# $par
# [1] 0.45118 0.24754 0.00000


