###################################################
### code chunk number 7: Cs1_Exercise1
###################################################
par(mfrow = c(3, 3))
sim.u <- -0.05
sim.Q <- 0.02
sim.R <- 0.05
nYr <- 50
fracmiss <- 0.1
init <- 7
years <- seq(1:nYr)
for (i in 1:9) {
  x <- rep(NA, nYr) # vector for ts w/o measurement error
  y <- rep(NA, nYr) # vector for ts w/ measurement error
  x[1] <- init
  for (t in 2:nYr) {
    x[t] <- x[t - 1] + sim.u + rnorm(1, mean = 0, sd = sqrt(sim.Q))
  }
  for (t in 1:nYr) {
    y[t] <- x[t] + rnorm(1, mean = 0, sd = sqrt(sim.R))
  }
  missYears <-
    sample(years[2:(nYr - 1)], floor(fracmiss * nYr), replace = FALSE)
  y[missYears] <- NA
  plot(years, y,
    xlab = "", ylab = "Log abundance", lwd = 2, bty = "l"
  )
  lines(years, x, type = "l", lwd = 2, lty = 2)
  title(paste("simulation ", i))
}
legend("topright", c("Observed", "True"),
  lty = c(-1, 2), pch = c(1, -1)
)


