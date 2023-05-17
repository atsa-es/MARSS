###################################################
### code chunk number 3: Cs02_plotwolfmoosedata
###################################################
x <- isleRoyal[, "Year"]
y <- log(isleRoyal[, c("Wolf", "Moose")])
graphics::matplot(x, y,
  ylab = "Log count", xlab = "Year", type = "l",
  lwd = 3, bty = "L", col = "black"
)
legend("topright", c("Wolf", "Moose"), lty = c(1, 2), bty = "n")


