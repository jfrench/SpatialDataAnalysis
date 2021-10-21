# understanding cepp.test internals
library(smerc)
# define information needed in function
data(nydf)
coords = as.matrix(nydf[,c("x", "y")])
cases = nydf$cases
pop = nydf$population
# expected number of cases in each region
ex = sum(cases)/sum(pop) * pop
# population upper bound
ubpop = 0.1
# don't use great circle distance
longlat = FALSE
# number of simulated data sets
nsim = 49
# simulation approach
simdist = "multinomial"
# statistic type
type = "poisson"

# number of observations
N = nrow(coords)
# distance between centroids
d = sp::spDists(coords, longlat = longlat)
# nn subject to population constraints
nn = scan.nn(d, pop, ubpop)
# cases in each candidate zone
yin = nn.cumsum(nn, cases)
# total number of cases
ty = sum(cases)
# expected cases in each candidate zone
ein = nn.cumsum(nn, ex)
# expected cases outside each candidate zone
eout = sum(ex) - ein
# observed test statistics for each candidate zone
tobs = stat.poisson(yin, ty - yin, ein, eout)

# max test statistic for simulated data sets
tsim = scan.sim(nsim = nsim, nn = nn, ty = ty, ex = ex,
                type = type, ein = ein, eout = eout, simdist = simdist,
                pop = pop)
# compute p values for each candidate zone
pvalue = mc.pvalue(tobs, tsim)

# determine most likely cluster information
mlc = which.max(tobs)
tobs[mlc]
pvalue[mlc]

# additional work has to be done to link the significant statistics
# with specific candidate zones
zone = nn2zones(nn)
zone[[mlc]]