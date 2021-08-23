# vector of distances between s0 and s1,...,s5
d = c(1, sqrt(2), 1, sqrt(1.25), 1)
# vector of responses
z = c(2, -3, 3, -4, 2)
# calculate weights
(w = d^(-2)/sum(d^(-2)))
# interpolated response
sum(z*w)
