# load packages
library(smacpod)

# load grave data
data(grave)

# perform test using spatial scan statistic
# use group name to select case group
scan = spscan.test(grave, nsim = 999, case = "affected")
# summary of scan test results
summary(scan)
# plot scan test results
plot(scan, chars = c(1, 20), main = "most likely cluster for grave data")
# extract most likely and other significant clusters
clusters(scan)

# q nearest neighbor test
# use position in levels(grave$marks) to select case group
qnn.test(grave, q = c(3, 5, 7, 9, 11, 13, 15), nsim = 499, case = 2)
