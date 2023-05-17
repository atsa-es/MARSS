###################################################
### code chunk number 4: Cs_030_plotdata
###################################################
par(mfrow = c(m, 1), mar = c(4, 4, 0.1, 0), oma = c(0, 0, 2, 0.5))
plot(years, dat, xlab = "", ylab = "Logit(s)", bty = "n", xaxt = "n", pch = 16, col = "darkgreen", type = "b")
plot(years, CUI.z, xlab = "", ylab = "CUI", bty = "n", xaxt = "n", pch = 16, col = "blue", type = "b")
axis(1, at = seq(1965, 2005, 5))
mtext("Year of ocean entry", 1, line = 3)


