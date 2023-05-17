###################################################
### code chunk number 4: fakepopdata
###################################################
x[1] <- init
for (t in 2:nYr) {
  x[t] <- x[t - 1] + sim.u + rnorm(1, mean = 0, sd = sqrt(sim.Q))
}


