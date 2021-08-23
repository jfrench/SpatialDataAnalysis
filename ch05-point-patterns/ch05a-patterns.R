library(splancs)
library(spatial)

set.seed(1) # for reproducability

# domain is unit square
domain = cbind(c(0,0,1,1), c(0,1,1,0))

# plots of clustering and regularity
# pdf(file = "csr.pdf", width = 7.5, height = 5)
par(mfrow = c(2,3), pty="s") #2x3 grid of plots
for (i in 1:6) {
  locs = csr(domain, 100) # generate date under CSR
  # plot points
  plot(locs[,1], locs[,2], xlab = "u", ylab = "v", cex.lab = 1.5, cex.axis = 1.1) 
  title("CSR")
}
# dev.off()

# pdf(file = "cluster_regular.pdf", width = 7.5, height = 5)
par (mfrow = c(2,3), pty="s") #2x3 grid of plots
for (i in 1:3) {
  test <- pcp.sim(rho = 5, m = 10, s2 = 0.01, region.poly = domain)
  plot(test, xlim = 0:1, ylim = 0:1, xlab = "u", ylab = "v", 
       cex.lab = 1.5, cex.axis = 1.1)
  title("Clustered")
}

# set boundaries of process to use Strauss function in spatial package
ppregion(xl = 0, xu = 1, yl = 0, yu = 1)
for (i in 1:3) {
  # Simulate 30 events from a Strauss process with probability c of other events
  # within distance r of other events.
  test <- Strauss(30, c = 0, r = 0.05)
  plot(test, xlim = 0:1, ylim = 0:1, xlab = "u", ylab = "v", cex.lab = 1.5, cex.axis = 1.1)
  title("Regular")
}
# dev.off()

clus = pcp.sim(rho = 1.4, m = 15, s2 = 0.01, region.poly = domain)
pp1 = vector("list", nrow(clus))
w = .03
for (i in 1:nrow(clus)) {
  ppregion(xl = clus[i,1] - w, 
           xu = clus[i,1] + w, 
           yl = clus[i,2] - w, 
           yu = clus[i,2] + w)
  temp = Strauss(5, c = 0, r = 0.027)
  pp1[[i]] = temp
}

# pdf("cluster_regular_events.pdf")
plot(0:1, 0:1, type = "n", xlab = "u", ylab = "v")
for (i in 1:nrow(clus)) points(pp1[[i]])
title("Cluster of Regular Event Locations")
# dev.off()

# pdf(file = "cluster_regular.pdf", width = 7.5, height = 5)
par(mfrow = c(2,3), pty = "s") #2x3 grid of plots
for (i in 1:3) {
  test <- pcp.sim(rho = 5, m = 10, s2 = 0.01, region.poly = domain)
  plot(test, xlim = 0:1, ylim = 0:1, 
       xlab = "u", ylab = "v", 
       cex.lab = 1.5, cex.axis = 1.1)
  title("Clustered")
}

# set boundaries of process to use Strauss function in spatial package
ppregion(xl = 0, xu = 1, yl = 0, yu = 1)
for (i in 1:3) {
  # Simulate 30 events from a Strauss process with probability c of other events
  # within distance r of other events.
  test <- Strauss(30, c = 0, r = 0.05)
  plot(test, xlim = 0:1, ylim = 0:1, 
       xlab = "u", ylab = "v", 
       cex.lab = 1.5, cex.axis = 1.1)
  title("Regular")
}
dev.off()

clus = pcp.sim(rho = 1.4, m = 5, s2 = 0.03, region.poly = domain)
pp1 = vector("list", nrow(clus))
w = .03
for (i in 1:nrow(clus)) {
  ppregion(xl = clus[i,1] - w, xu = clus[i,1] + w, 
           yl = clus[i,2] - w, yu = clus[i,2] + w)
  temp = Strauss(5, c = 0, r = 0.027)
  pp1[[i]] = temp
}

# pdf("cluster_regular_events.pdf")
plot(0:1, 0:1, type = "n", xlab = "u", ylab = "v")
for (i in 1:nrow(clus)) points(pp1[[i]])
title("Cluster of Regular Event Locations")
# dev.off()

ppregion(xl = 0, xu = 1, yl = 0, yu = 1)
reg = Strauss(5, c = 0, r = .5)
reg = cbind(reg$x, reg$y)
pp2 = vector("list", nrow(reg))
w = .05
for (i in 1:nrow(reg)) {
  pt = reg[i,]
  dom = cbind(c(pt[1] - w, pt[1] - w, pt[1] + w, pt[1] + w), 
              c(pt[2] - w, pt[2] + w, pt[2] + w, pt[2] - w))
  pp2[[i]] = csr(dom, 10)
}
# pdf("regular_cluster_events.pdf")
plot(0:1, 0:1, type = "n", xlab = "u", ylab = "v")
for (i in 1:nrow(reg)) points(pp2[[i]])
title("Regular cluster of event locations")
# dev.off()

# figures 5.5 and 5.6
x <- 1:20
y <- 1:20

z <- matrix(0,20,20)
mu1 <- c(3,3)
mu2 <- c(16,14)
sigmasq1 <- 6
sigmasq2 <- 12
Siginv <- matrix(0,2,2)
Siginv[1,1] <- 1/sigmasq1
Siginv[2,2] <- 1/sigmasq2

Siginv2 <- matrix(0,2,2)
Siginv2[1,1] <- 1/60
Siginv2[2,2] <- 1/35

for (i in 1:20) {
  for (j in 1:20) {
    z[i,j] <- (1/(2*pi))*exp((-1/2)*
                               ( (c(x[i],y[j]) - mu1)%*%Siginv%*%(c(x[i],y[j]) - mu1) ) )+
      (1/(2*pi))*exp((-1/2)*
                       ( (c(x[i],y[j]) - mu2)%*%Siginv2%*%(c(x[i],y[j]) - mu2) ) )
  }
  print(i)
}

par(mfrow=c(1,2),pty="s")
persp(x,y,z,theta=45,phi=45,xlab="u",ylab="v",zlab="lambda",cex.lab=1.5)
contour(x,y,z,xlim=c(0,20),ylim=c(0,20),xlab="u",ylab="v",cex.lab=1.5,labcex=1.2,cex=1.5,cex.axis=1.5)

###################################
# Generate 200 events following a
# homogeneous Poisson process (CSR).
# We will thin these down to 100
# events following the heterogeneous
# process.
################################

xrand <- runif(200,min=0,max=20)
yrand <- runif(200,min=0,max=20)

test <- runif(200,min=,max=max(z))

# Define the z-value at each point corresponding to the
# defined intensity function.
zval <- 1:200
for (i in 1:200) {
  zval[i] <- z[ trunc(xrand[i])+1,trunc(yrand[i])+1 ]
}

####################################
# Make six examples (using a rejection sampling method).  Keep the first 100 points with
# associated z-value below the density value.
####################################

xrand <- runif(200,min=0,max=20)
yrand <- runif(200,min=0,max=20)
test <- runif(200,min=,max=max(z))
zval <- 1:200
for (i in 1:200) {
  zval[i] <- z[trunc(xrand[i])+1,trunc(yrand[i])+1]
}
plot1x <- xrand[test<zval][1:100]
plot1y <- yrand[test<zval][1:100]

xrand <- runif(200,min=0,max=20)
yrand <- runif(200,min=0,max=20)
test <- runif(200,min=,max=max(z))
zval <- 1:200
for (i in 1:200) {
  zval[i] <- z[ trunc(xrand[i])+1,trunc(yrand[i])+1 ]
}
plot2x <- xrand[test<zval][1:100]
plot2y <- yrand[test<zval][1:100]

xrand <- runif(200,min=0,max=20)
yrand <- runif(200,min=0,max=20)
test <- runif(200,min=,max=max(z))
zval <- 1:200
for (i in 1:200) {
  zval[i] <- z[ trunc(xrand[i])+1,trunc(yrand[i])+1 ]
}
plot3x <- xrand[test<zval][1:100]
plot3y <- yrand[test<zval][1:100]

xrand <- runif(200,min=0,max=20)
yrand <- runif(200,min=0,max=20)
test <- runif(200,min=,max=max(z))
zval <- 1:200
for (i in 1:200) {
  zval[i] <- z[ trunc(xrand[i])+1,trunc(yrand[i])+1 ]
}
plot4x <- xrand[test<zval][1:100]
plot4y <- yrand[test<zval][1:100]

xrand <- runif(200,min=0,max=20)
yrand <- runif(200,min=0,max=20)
test <- runif(200,min=,max=max(z))
zval <- 1:200
for (i in 1:200) {
  zval[i] <- z[ trunc(xrand[i])+1,trunc(yrand[i])+1 ]
}
plot5x <- xrand[test<zval][1:100]
plot5y <- yrand[test<zval][1:100]

xrand <- runif(200,min=0,max=20)
yrand <- runif(200,min=0,max=20)
test <- runif(200,min=,max=max(z))
zval <- 1:200
for (i in 1:200) {
  zval[i] <- z[ trunc(xrand[i])+1,trunc(yrand[i])+1 ]
}
plot6x <- xrand[test<zval][1:100]
plot6y <- yrand[test<zval][1:100]

par(mfrow=c(2,3))

plot(plot1x,plot1y,xlim=c(0,20),ylim=c(0,20),xlab="u",ylab="v")
contour(x,y,z,lwd=2,lty=3,add=T,drawlabels=FALSE)

plot(plot2x,plot2y,xlim=c(0,20),ylim=c(0,20),xlab="u",ylab="v")
contour(x,y,z,lwd=2,lty=3,add=T,drawlabels=FALSE)

plot(plot3x,plot3y,xlim=c(0,20),ylim=c(0,20),xlab="u",ylab="v")
contour(x,y,z,lwd=2,lty=3,add=T,drawlabels=FALSE)

plot(plot4x,plot4y,xlim=c(0,20),ylim=c(0,20),xlab="u",ylab="v")
contour(x,y,z,lwd=2,lty=3,add=T,drawlabels=FALSE)

plot(plot5x,plot5y,xlim=c(0,20),ylim=c(0,20),xlab="u",ylab="v")
contour(x,y,z,lwd=2,lty=3,add=T,drawlabels=FALSE)

plot(plot6x,plot6y,xlim=c(0,20),ylim=c(0,20),xlab="u",ylab="v")
contour(x,y,z,lwd=2,lty=3,add=T,drawlabels=FALSE)
