###################################################
### code chunk number 37: Cs_410_marssperffig
###################################################
if (sims.exist) { # flag to check that the code to creat the plot has been run
  # This makes a plot of the comparisons
  par(ylbias = -.2, tcl = -.2, cex = .75)
  graphics::boxplot(t(params[c(7, 1, 3, 5, 8, 2, 4, 6), ]), names = c("AR(2)\n", "MARSS\nEM", "MARSS\nBFGS", "ARMA\n(2,2)", "AR(2)\n", "MARSS\nEM", "MARSS\nBFGS", "ARMA\n(2,2)"), ylab = "estimates of the ar coefficients", las = 2)
  points(1:8, apply(params[c(7, 1, 3, 5, 8, 2, 4, 6), ], 1, mean), pch = "x", cex = 1.25)
  par(cex = 1.5)
  axis(side = 3, at = c(2, 6), labels = c(expression(b[1]), expression(b[2])), tick = FALSE, cex = 2)
  lines(c(0, 4.5), c(true.2ss[2], true.2ss[2]), lty = 2)
  lines(c(4.5, 9), c(true.2ss[3], true.2ss[3]), lty = 2)
  abline(v = 4.5)
}


