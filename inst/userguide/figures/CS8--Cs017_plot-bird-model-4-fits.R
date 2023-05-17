###################################################
### code chunk number 33: Cs017_plot-bird-model-4-fits
###################################################
# Make a plot of the predicted trajectory, confidence intervals, and the raw data in log-space
plot(kestrel[, 1], kestrel[, 2], xlab = "", ylab = "Index of kestrel abundance", 
     main = "", col = "red", ylim = c(0, 2), pch = 21)
points(kestrel[, 1], kestrel[, 3], col = "blue", pch = 22)
points(kestrel[, 1], kestrel[, 4], col = "purple", pch = 25)
lines(kestrel[, 1], c(kem.b4$states[1, ]), lty = 3, lwd = 2, col = "red")
lines(kestrel[, 1], c(kem.b4$states[2, ]), lty = 3, lwd = 2, col = "blue")
lines(kestrel[, 1], c(kem.b4$states[2, ] + coef(kem.b4, type = "matrix")$A[3, 1]), 
      lty = 3, lwd = 2, col = "purple")
legend("topright", inset = 0.1, legend = c("British Columbia", "Alberta", "Saskatchewan"), 
       col = c("red", "blue", "purple"), pch = c(21, 22, 25))


