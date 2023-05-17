###################################################
### code chunk number 12: Cs_10_plot-dlm-Forecast-Logit
###################################################
par(mar = c(4, 4, 0.1, 0), oma = c(0, 0, 2, 0.5))
ylims <- c(min(fore.mean - 2 * sqrt(fore.var)), max(fore.mean + 2 * sqrt(fore.var)))
plot(years, t(dat),
  type = "p", pch = 16, ylim = ylims,
  col = "blue", xlab = "", ylab = "Logit(s)", xaxt = "n"
)
lines(years, fore.mean, type = "l", xaxt = "n", ylab = "", lwd = 3)
lines(years, fore.mean + 2 * sqrt(fore.var))
lines(years, fore.mean - 2 * sqrt(fore.var))
axis(1, at = seq(1965, 2005, 5))
mtext("Year of ocean entry", 1, line = 3)


