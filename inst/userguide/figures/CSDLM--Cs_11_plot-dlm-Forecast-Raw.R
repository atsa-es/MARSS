###################################################
### code chunk number 13: Cs_11_plot-dlm-Forecast-Raw
###################################################
invLogit <- function(x) {
  1 / (1 + exp(-x))
}
ff <- invLogit(fore.mean)
fup <- invLogit(fore.mean + 2 * sqrt(fore.var))
flo <- invLogit(fore.mean - 2 * sqrt(fore.var))
par(mar = c(4, 4, 0.1, 0), oma = c(0, 0, 2, 0.5))
ylims <- c(min(flo), max(fup))
plot(years, invLogit(t(dat)),
  type = "p", pch = 16, ylim = ylims,
  col = "blue", xlab = "", ylab = "Survival", xaxt = "n"
)
lines(years, ff, type = "l", xaxt = "n", ylab = "", lwd = 3)
lines(years, fup)
lines(years, flo)
axis(1, at = seq(1965, 2005, 5))
mtext("Year of ocean entry", 1, line = 3)


