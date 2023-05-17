###################################################
### code chunk number 28: Cs_402_fig-arss-model
###################################################
TT <- 1000 # set long
true.2ss <- c(r = .5, b1 = -1.5, b2 = -0.75, q = .1)
temp <- arima.sim(
  n = TT, list(ar = true.2ss[c("b1", "b2")]),
  sd = sqrt(true.2ss["q"])
)
sim.ar <- matrix(temp, nrow = 1)
noise <- rnorm(TT - 1, 0, sqrt(true.2ss["r"]))
noisy.data <- sim.ar[2:TT] + noise


