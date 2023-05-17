###################################################
### code chunk number 27: Cs008_plot-bird-data
###################################################
# Make a plot of the three time series
plot(kestrel[, 1], kestrel[, 2], xlab = "", ylab = "Index of kestrel abundance", main = "", col = "red", ylim = c(0, 2), pch = 21)
points(kestrel[, 1], kestrel[, 3], col = "blue", pch = 22)
points(kestrel[, 1], kestrel[, 4], col = "purple", pch = 25)
legend("topright", inset = 0.1, legend = c("British Columbia", "Alberta", "Saskatchewan"), col = c("red", "blue", "purple"), pch = c(21, 22, 25))


