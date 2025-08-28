# N = 10
# test statistic is (sample mean - 0)/(sd(data)/sqrt(N))

set.seed(1)

Tobs = 4
N = 10 # sample size
nsim = 999 # number of simulated data sets

# Under H0, say the distribution of the data is N(0, 1)
# Simulate data sets under H0
simdata = matrix(rnorm(N * nsim, mean = 0, sd = 1),
                 nrow = nsim, ncol = N)
# calculate test statistic for each simulated data set
Tsim = rowMeans(simdata)/(apply(simdata, 1, sd)/sqrt(N))

# plot approximate null distribution
plot(density(Tsim))
abline(v = 4) # location of Tobs

# include observed test statistic to simulated test statistics
Tsim = c(Tsim, Tobs)

# determine proportion of test statistics as extreme as the observed
# test statistic
mean(Tsim >= Tobs)
