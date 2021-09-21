# install.packages(c("SpatialTools", "KernSmooth")
library(smerc)

############################################
# Apply Tango's index to NY leukemia data
############################################
# read in ny data
nydf <- read.table("./data/NYTRACT.dat")
names(nydf) = c("x", "y", "pop", "cases")

###################################################
coords = as.matrix(nydf[,c("x", "y")])
cases = nydf$cases
pop = nydf$pop

# Find distance matrix
d = as.matrix(dist(coords))

##################################################
# Exponential decay weight matrix
# use different kappas in defining weights
w1r  <- dweights(coords, kappa = 1, pop = pop, type = "rogerson")
w7r  <- dweights(coords, kappa = 7, pop = pop, type = "rogerson")
w10r <- dweights(coords, kappa = 10, pop = pop, type = "rogerson")
w15r <- dweights(coords, kappa = 15, pop = pop, type = "rogerson")
w20r <- dweights(coords, kappa = 20, pop = pop, type = "rogerson")

# Code to plot weights and "effective range"
# plotdist <- 0:max(d)
# kappa <- 10
# plot(plotdist,
#      exp(-plotdist/kappa),type="l",xlab="Distance",ylab="exp(-distance/kappa)",
#     cex.lab=1.5,cex.axis=1.25,ylim=c(0,1))
# #rug(dist)
# title(paste("kappa = ",kappa),cex.main=2.0)
# effrange <- -kappa*log(0.05)
# segments(0,0.05,effrange,0.05)
# segments(effrange,0,effrange,0.05)

###################################
# Calculate Tango's statistic

# the tango function takes the number of cases,
# the population, and the
# matrix of weights
(tango_1  <- tango.test(cases, pop, w1r))
(tango_7  <- tango.test(cases, pop, w7r))
(tango_10 <- tango.test(cases, pop, w10r))
(tango_15 <- tango.test(cases, pop, w15r))
(tango_20 <- tango.test(cases, pop, w20r))

# extracting goodness-of-fit and spatial autocorrelation
# components of tango's index
gof <- c(tango_1$gof,tango_7$gof,
         tango_10$gof,tango_15$gof,
         tango_20$gof)
sa <- c(tango_1$sa,tango_7$sa,
        tango_10$sa,tango_15$sa,
        tango_20$sa)
plot(gof, sa)
# gof stays the same since all sets of weights have w_{ii} =
# 1 the difference between value's of tango's index for
# different values of kappa derives entirely from
# differences in the spatial autocorrelation component
#
# changing wstar changes the skewness of tstat, which reduces
# the associated df for the chi-square approximation,
# resulting in much smaller p-values

###################################
# Monte Carlo p-values

# compare monte carlo p-value to chi-square approximation p-value
(tango_mc1 <-  tango.test(cases, pop, w1r, nsim = 9999))
(tango_mc7 <-  tango.test(cases, pop, w7r, nsim = 9999))
(tango_mc10 <-  tango.test(cases, pop, w10r, nsim = 9999))
(tango_mc15 <-  tango.test(cases, pop, w15r, nsim = 9999))
(tango_mc20 <-  tango.test(cases, pop, w20r, nsim = 9999))

# comparing gof and sa components of tango's index for the observed
# data to the simulated data
# x is observed
plot(tango_mc1)
plot(tango_mc7)
plot(tango_mc10)
plot(tango_mc15)
plot(tango_mc20)
