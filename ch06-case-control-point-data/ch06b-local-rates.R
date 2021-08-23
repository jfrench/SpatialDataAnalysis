# load packages
library(smacpod)

# load grave data
data(grave)

# perform test using spatial scan statistic
scan = spscan.test(grave, nsim = 999, case = 2)
plot(scan, chars = c(1, 20), main = "most likely cluster for grave data")

# q nearest neighbor test
qnn.test(grave, q = c(3, 5, 7, 9, 11, 13, 15), nsim = 499, case = 2)
