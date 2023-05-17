###################################################
### code chunk number 22: Covar_sec6_12_plot-seas-effects
###################################################
par(mfrow = c(3, 1), mar = c(2, 4, 2, 2))
graphics::matplot(t(seas.1), type = "l", bty = "n", xaxt = "n", ylab = "Fixed monthly", col = 1:5)
axis(1, labels = month.abb, at = 1:12, las = 2, cex.axis = 0.75)
legend("topright", lty = 1:5, legend = phytos, cex = 0.6, col = 1:5)

graphics::matplot(t(seas.2), type = "l", bty = "n", xaxt = "n", ylab = "Cubic", col = 1:5)
axis(1, labels = month.abb, at = 1:12, las = 2, cex.axis = 0.75)
legend("topright", lty = 1:5, legend = phytos, cex = 0.6, col = 1:5)

graphics::matplot(t(seas.3), type = "l", bty = "n", xaxt = "n", ylab = "Fourier", col = 1:5)
axis(1, labels = month.abb, at = 1:12, las = 2, cex.axis = 0.75)
legend("topright", lty = 1:5, legend = phytos, cex = 0.6, col = 1:5)


