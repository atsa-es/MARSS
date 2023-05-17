###################################################
### code chunk number 19: plotslopetests-www
###################################################
library(Hmisc)
resids <- MARSSresiduals(kem.3, type = "tT")$mar.residuals
x <- seq(tsp(WWWusage)[1], tsp(WWWusage)[2], tsp(WWWusage)[3])
plot(x, resids[3, ],
  ylab = "Mar. residuals", xlab = "", type = "l", ylim = c(-4, 4),
  xaxp = c(tsp(WWWusage)[1] - 1, tsp(WWWusage)[2], 10)
)
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
abline(h = c(1.97, -1.97), lty = 2)
abline(h = 0)
title("test for slope changes")


