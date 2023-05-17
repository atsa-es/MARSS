###################################################
### code chunk number 23: Cs_301_sim-ar3-data
###################################################
TT <- 100
true.3 <- c(r = 0, b1 = -1.5, b2 = -0.75, b3 = .05, q = 1)
temp3 <- arima.sim(
  n = TT, list(ar = true.3[c("b1", "b2", "b3")]),
  sd = sqrt(true.3["q"])
)
sim.ar3 <- matrix(temp3, nrow = 1)


