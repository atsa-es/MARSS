###################################################
### code chunk number 14: Cs11_plotresids
###################################################
par(mfrow = c(2, 1), mar = c(4, 3, 2, 1))
x <- seq(tsp(Nile)[1], tsp(Nile)[2], tsp(Nile)[3])
plot(x, resids.2[1, ], ylab = "", xlab = "", type = "l", ylim = c(-4, 4), xaxp = c(1870, 1970, 10))
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
abline(h = c(1.97, -1.97), lty = 2)
title("test for outliers")

plot(x, resids.2[2, ], ylab = "", xlab = "", type = "l", ylim = c(-4, 4), xaxp = c(1870, 1970, 10))
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
abline(h = c(1.97, -1.97), lty = 2)
title("test for level changes")
mtext("Standardized residuals", side = 2, outer = TRUE, line = -1)


