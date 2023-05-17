###################################################
### code chunk number 36: Cs708_multivariate
###################################################
par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
covariate <- matrix(c(rep(0, TT - 10), rep(1, 10)), nrow = 1)
ymc <- ym
D <- matrix(c(-1, -1, 0, 1, 1), ncol = 1)
ymc <- ym + D %*% covariate
matplot(t, t(ymc), pch = 1:n, ylab = "y", xlab = "", main = "data")
lines(level)
plot(t, covariate[1, ], col = "blue", lty = 2, type = "l", main = "covariate")


