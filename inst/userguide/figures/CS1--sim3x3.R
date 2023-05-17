###################################################
### code chunk number 6: sim3x3
###################################################
par(mfrow = c(3, 3))
sim.u <- -0.05 # growth rate
sim.Q <- 0.02 # process error variance
sim.R <- 0.05 # non-process error variance
nYr <- 50 # number of years of data to generate
fracmissing <- 0.1 # fraction of years that are missing
init <- 7 # log of initial pop abundance

years <- seq(1:nYr) # col of years
x.tss <- matrix(NA, nrow = 9, ncol = nYr) # creates vector for ts w/o obs. error
y.tss <- matrix(NA, nrow = 9, ncol = nYr) # creates vector for ts w/ obs. error
for (i in 1:9) {
  x.tss[i, 1] <- init
  for (t in 2:nYr) x.tss[i, t] <- x.tss[i, t - 1] + sim.u + rnorm(1, mean = 0, sd = sqrt(sim.Q))
  for (t in 1:nYr) y.tss[i, t] <- x.tss[i, t] + rnorm(1, mean = 0, sd = sqrt(sim.R))
  missYears <- sample(years[2:(nYr - 1)], floor(fracmissing * nYr), replace = FALSE)
  y.tss[i, missYears] <- NA
  plot(years, y.tss[i, ], xlab = "", ylab = "Index of log abundance", lwd = 2, bty = "l")
  lines(years, x.tss[i, ], type = "l", lwd = 2, lty = 2)
  title(paste("simulation ", i))
}
# legend("topright", c("Observed","True"),lty = c(-1, 2), pch = c(1, -1))


