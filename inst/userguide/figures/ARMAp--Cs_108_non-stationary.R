###################################################
### code chunk number 11: Cs_108_non-stationary
###################################################
TT <- 50
true.2 <- c(r = 0, b1 = -1.5, b2 = -0.75, q = 1)
sim.ar2.ns <- rep(NA, TT)
sim.ar2.ns[1] <- -30
sim.ar2.ns[2] <- -10
for (i in 3:TT) {
  sim.ar2.ns[i] <- true.2[2] * sim.ar2.ns[i - 1] +
    true.2[3] * sim.ar2.ns[i - 2] + rnorm(1, 0, sqrt(true.2[4]))
}

model.list.3 <- model.list.2
model.list.3$x0 <- matrix(sim.ar2.ns[2:1], 2, 1)

ar3.marss <- MARSS(sim.ar2.ns[3:TT], model = model.list.3, silent = TRUE)
ar3.arima <- arima(sim.ar2.ns[3:TT], order = c(2, 0, 0), include.mean = FALSE)

print(cbind(
  true = true.2[2:4],
  estimates.marss = coef(ar3.marss, type = "vector"),
  estimates.arima = c(coef(ar3.arima, type = "vector"), ar3.arima$sigma2)
))


