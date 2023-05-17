###################################################
### code chunk number 16: Cs13_plot-plankton-data
###################################################
graphics::matplot((1:(52 * 6))[27:295], t(d.plank.dat), type = "l", lty = c(1, 1, 1, 1), lwd = c(1, 1, 3, 3), xlab = "week of study", ylab = "log biomass", xaxt = "n", xlim = c(11, 52 * 6 - 11), bty = "L")
# axis(1,at=(1:(52*6))[seq(27,295,2)])
axis(1, at = seq(1, 52 * 6, 2))
abline(v = c(52 * (1:6)))
abline(h = 0)


