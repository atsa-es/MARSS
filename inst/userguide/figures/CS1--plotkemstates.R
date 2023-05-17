###################################################
### code chunk number 11: plotkemstates
###################################################
par(mfrow = c(3, 3))
for (i in 1:9) {
  plot(years, y.tss[i, ], xlab = "", ylab = "Index of log abundance", lwd = 2, bty = "l")
  lines(years, x.tss[i, ], type = "l", lwd = 2, lty = 2)
  lines(years, kem.states[i, ], type = "l", col = 2, lwd = 1, lty = 1)
  title(paste("simulation ", i))
}
# legend("topright", c("Observed","True","KalmanEM estimate"),lty = c(-1, 2, 1), pch = c(1, -1, -1),col=c(1,1,2))


