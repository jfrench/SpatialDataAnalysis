h <- seq(0, 1, len = 1000)

vsph <- 1.25 - cov.spatial(h, cov.model = "spherical",
                    cov.pars = c(1, 0.75))

png("sph.png", width = 6, height = 5, units = "in", res = 300)
plot(h, vsph, type = "l",
     main = "spherical", ylab = "semivariance",
     ylim = c(0, 1.5))
dev.off()

vexp <- 1.25 - cov.spatial(h, cov.model = "exponential",
                           cov.pars = c(1, 0.25))

png("exp.png", width = 6, height = 5, units = "in", res = 300)
plot(h, vexp, type = "l",
     main = "exponential", ylab = "semivariance",
     ylim = c(0, 1.5))
dev.off()

vgau <- 1.25 - cov.spatial(h, cov.model = "gaussian",
                           cov.pars = c(1, 0.75/sqrt(3)))

png("gau.png", width = 6, height = 5, units = "in", res = 300)
plot(h, vgau, type = "l", ylab = "semivariance",
     main = "gaussian",
     ylim = c(0, 1.5))
dev.off()

vmata <- 1.25 - cov.spatial(h, cov.model = "matern",
                            cov.pars = c(1, 0.25),
                            kappa = 0.5)

vmatb <- 1.25 - cov.spatial(h, cov.model = "matern",
                            cov.pars = c(1, 0.25),
                            kappa = 1.5)

vmatc <- 1.25 - cov.spatial(h, cov.model = "matern",
                            cov.pars = c(1, 0.25),
                            kappa = 2.5)

vmat <- cbind(vmata, vmatb, vmatc)

png("matern.png", width = 6, height = 5, units = "in", res = 300)
matplot(h, vmat, type = "l", ylab = "semivariance",
        col = c("black", "red", "blue"))
legend("bottomright",
       legend = c(expression(alpha == 0.5),
                  expression(alpha == 1.5),
                  expression(alpha == 2.5)),
       lty = 1:3, col = c("black", "red", "blue"))
dev.off()
