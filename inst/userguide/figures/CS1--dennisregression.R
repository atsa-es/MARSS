###################################################
### code chunk number 14: dennisregression
###################################################
par(mfrow = c(1, 1))
plot(tau, delta.pop, xlab = "Time step size (tau)", ylab = "Population transition size", xlim = c(0, max(tau)), bty = "l")
abline(den91, col = 2)


