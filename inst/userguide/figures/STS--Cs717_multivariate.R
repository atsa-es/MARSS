###################################################
### code chunk number 45: Cs717_multivariate
###################################################
par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
ylims <- c(min(ym, na.rm = TRUE), max(ym, na.rm = TRUE))
plot(t, ym[1, ], ylim = ylims, type = "p")
lines(t, level1, col = "black")
plot(t, ym[2, ], ylim = ylims, type = "p")
lines(t, level2, col = "black")


