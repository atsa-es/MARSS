###################################################
### code chunk number 9: Cs_07_plotdlm1
###################################################
ylabs <- c(expression(alpha[t]), expression(beta[t]))
colr <- c("darkgreen", "blue")
par(mfrow = c(m, 1), mar = c(4, 4, 0.1, 0), oma = c(0, 0, 2, 0.5))
for (i in 1:m) {
  mn <- dlm1$states[i, ]
  se <- dlm1$states.se[i, ]
  plot(years, mn,
    xlab = "", ylab = ylabs[i], bty = "n", xaxt = "n", type = "n",
    ylim = c(min(mn - 2 * se), max(mn + 2 * se))
  )
  lines(years, rep(0, TT), lty = "dashed")
  lines(years, mn, col = colr[i], lwd = 3)
  lines(years, mn + 2 * se, col = colr[i])
  lines(years, mn - 2 * se, col = colr[i])
}
axis(1, at = seq(1965, 2005, 5))
mtext("Year of ocean entry", 1, line = 3)


