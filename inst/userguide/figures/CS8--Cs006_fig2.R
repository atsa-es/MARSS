###################################################
### code chunk number 9: Cs006_fig2
###################################################
# Code for plotting the fit from the best model
plot(okanaganRedds[, 1], logRedds[1, ], xlab = "", ylab = "Redd counts", 
     main = "", col = "red", ylim = c(0, 8)
)
points(okanaganRedds[, 1], logRedds[2, ], col = "blue", pch = 2)
lines(okanaganRedds[, 1], c(kem1$states), lty = 1, lwd = 2)
lines(okanaganRedds[, 1], c(kem1$states + 2 * kem1$states.se), lty = 1, lwd = 1, col = "grey40")
lines(okanaganRedds[, 1], c(kem1$states - 2 * kem1$states.se), lty = 1, lwd = 1, col = "grey40")


