###################################################
### code chunk number 30: Cs702_multivariate
###################################################
par(mfrow = c(2, 1), mar = c(3, 3, 1, 1))
ylims <- c(min(ym, na.rm = TRUE), max(ym, na.rm = TRUE))
plot(t, trend, ylim = ylims, col = "red", type = "l")
lines(t, level, col = "black")
legend("topright", c("trend", "level"), lty = 1, col = c("red", "black"))
matplot(t, t(ym), pch = 1:n, ylab = "y", xlab = "", ylim = ylims, main = "bad data")
lines(t, level)


