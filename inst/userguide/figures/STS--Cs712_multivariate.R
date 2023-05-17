###################################################
### code chunk number 40: Cs712_multivariate
###################################################
r2 <- r * c(100, 10, 10, 200, 400)
a <- runif(n, -1, 1)
err <- rnorm(n * TT, mean = rep(a, each = TT), sd = rep(sqrt(r2), each = TT))
ym2 <- matrix(1, nrow = n) %*% level + matrix(err, nrow = n, byrow = TRUE)
ym2[sample(n * TT, miss.percent * n * TT)] <- NA
matplot(t, t(ym2), pch = 1:n, ylab = "y", xlab = "", main = "data with different error and bias")
lines(level)


