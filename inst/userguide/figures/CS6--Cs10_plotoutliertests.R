###################################################
### code chunk number 13: Cs10_plotoutliertests
###################################################
library(Hmisc)
par(mfrow = c(3, 1), mar = c(3, 4, 1.5, 2))
x <- seq(tsp(Nile)[1], tsp(Nile)[2], tsp(Nile)[3])
plot(x, resids.0[1, ],
  ylab = "Std. residuals", xlab = "", type = "l",
  ylim = c(-4, 4), xaxp = c(1870, 1970, 10), bty = "L"
)
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
abline(h = c(1.97, -1.97, 0), lty = 2)
title("model 0--flat level")

plot(x, resids.1[1, ],
  ylab = "Std. residuals", xlab = "", type = "l",
  ylim = c(-4, 4), xaxp = c(1870, 1970, 10), bty = "L"
)
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
abline(h = c(1.97, -1.97, 0), lty = 2)
title("model 1--linearly declining level")

plot(x, resids.2[1, ],
  ylab = "Std. residuals", xlab = "", type = "l",
  ylim = c(-4, 4), xaxp = c(1870, 1970, 10), bty = "L"
)
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
abline(h = c(1.97, -1.97, 0), lty = 2)
title("model 2--stochastic level")


