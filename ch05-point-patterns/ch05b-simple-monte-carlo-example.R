# N = 10
# test statistic is (sample mean - 0)/(sd(data)/sqrt(N))

Tobs = 4 
# Under H0, say the distribution of the data is N(0, 1)
# Simulate data sets under H0
simdata = matrix(rnorm(10 * 999, mean = 0, sd = 1), nrow = 999, ncol = 10)
# calculate test statistic for each simulated data set
Tsim = rowMeans(simdata)/(apply(simdata, 1, sd)/sqrt(10))
# include observed test statistic to simulated test statistics
Tsim = c(Tsim, Tobs)
# determine proportion of test statistics as extreme as the observed
# test statistc
mean(Tsim >= Tobs)
