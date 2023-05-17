###################################################
### code chunk number 11: Cs08_plotfit
###################################################
library(Hmisc)
par(mfrow = c(3, 1), mar = c(4, 4, 0.5, 0.5), oma = c(1, 1, 1, 1))
x <- seq(tsp(Nile)[1], tsp(Nile)[2], tsp(Nile)[3])
# model 0
plot(Nile, ylab = "Flow volume", xlab = "", xaxp = c(1870, 1970, 10), bty = "L")
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
kem <- kem.0 # model 0 results
lines(x, kem$states[1, ], col = "red", lwd = 2)
legend("topright", paste("model 0, AICc=", format(kem.0$AICc, digits = 1)), bty = "n")

# model 1
plot(Nile, ylab = "Flow volume", xlab = "", xaxp = c(1870, 1970, 10), bty = "n")
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
kem <- kem.1 # model 1 results
lines(x, kem$states[1, ], col = "red", lwd = 2)
legend("topright", paste("model 1, AICc=", format(kem.1$AICc, digits = 1)), bty = "n")

# model 2
plot(Nile, ylab = "Flow volume", xlab = "", xaxp = c(1870, 1970, 10), bty = "L")
minor.tick(nx = 10, ny = 0, tick.ratio = .3)
kem <- kem.2 # model 0 results
lines(x, kem$states[1, ], col = "red", lwd = 2)
lines(1871:1970, kem$states[1, ] - 2 * kem$states.se[1, ], col = "red", lty = 2)
lines(1871:1970, kem$states[1, ] + 2 * kem$states.se[1, ], col = "red", lty = 2)
legend("topright", paste("model 2, AICc=", format(kem$AICc, digits = 1)), bty = "n")


