###################################################
### code chunk number 4: Cs_101_ar2-sim
###################################################
TT <- 50
true.2 <- c(r = 0, b1 = -1.5, b2 = -0.75, q = 1)
temp <- arima.sim(n = TT, list(ar = true.2[2:3]), sd = sqrt(true.2[4]))
sim.ar2 <- matrix(temp, nrow = 1)


