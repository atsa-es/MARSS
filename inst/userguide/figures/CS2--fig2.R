###################################################
### code chunk number 9: fig2
###################################################
kem1 <- MARSS(dat,
  model =
    list(Z = factor(c(1, 1, 1, 1, 1)), R = "diagonal and equal")
)
graphics::matplot(years, t(dat), xlab = "", ylab = "Index of log abundance", pch = c("1", "2", "3", "4", "5"), ylim = c(5, 9), bty = "L")
lines(years, kem1$states - 1.96 * kem1$states.se, type = "l", lwd = 1, lty = 2, col = "red")
lines(years, kem1$states + 1.96 * kem1$states.se, type = "l", lwd = 1, lty = 2, col = "red")
lines(years, kem1$states, type = "l", lwd = 2)
title("Observations and total population estimate", cex.main = .9)


