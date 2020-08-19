###################################################
### code chunk number 2: Cs_000_required_libraries
###################################################
library(MARSS)


###################################################
### code chunk number 4: Cs_101_ar2-sim
###################################################
TT <- 50
true.2 <- c(r = 0, b1 = -1.5, b2 = -0.75, q = 1)
temp <- arima.sim(n = TT, list(ar = true.2[2:3]), sd = sqrt(true.2[4]))
sim.ar2 <- matrix(temp, nrow = 1)


###################################################
### code chunk number 5: Cs_102_ar2-model
###################################################
Z <- matrix(c(1, 0), 1, 2)
B <- matrix(list("b1", 1, "b2", 0), 2, 2)
U <- matrix(0, 2, 1)
Q <- matrix(list("q", 0, 0, 0), 2, 2)
A <- matrix(0, 1, 1)
R <- matrix(0, 1, 1)
mu <- matrix(sim.ar2[2:1], 2, 1)
V <- matrix(0, 2, 2)
model.list.2 <- list(Z = Z, B = B, U = U, Q = Q, A = A, 
                     R = R, x0 = mu, V0 = V, tinitx = 0)


###################################################
### code chunk number 6: Cs_103_ar2-fit
###################################################
ar2 <- MARSS(sim.ar2[3:TT], model = model.list.2)


###################################################
### code chunk number 7: Cs_104_ar2-fit
###################################################
print(cbind(true = true.2[2:4], estimates = coef(ar2, type = "vector")))


###################################################
### code chunk number 8: Cs_105_ar2-gappy
###################################################
gappy.data <- sim.ar2[3:TT]
gappy.data[floor(runif(TT / 2, 2, TT))] <- NA
ar2.gappy <- MARSS(gappy.data, model = model.list.2)


###################################################
### code chunk number 9: Cs_106_ar2-gappy
###################################################
print(cbind(
  true = true.2[2:4],
  estimates.no.miss = coef(ar2, type = "vector"),
  estimates.w.miss = coef(ar2.gappy, type = "vector")
))


###################################################
### code chunk number 10: Cs_107_arima
###################################################
arima(gappy.data, order = c(2, 0, 0), include.mean = FALSE)


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


###################################################
### code chunk number 13: Cs_201_mar2-sim
###################################################
TT <- 50
true.2 <- c(r = 0, b1 = -1.5, b2 = -0.75, q = 1)
temp1 <- arima.sim(n = TT, list(ar = true.2[c("b1", "b2")]), sd = sqrt(true.2["q"]))


###################################################
### code chunk number 15: Cs_202_mar2-sim
###################################################
temp2 <- arima.sim(n = TT, list(ar = true.2[c("b1", "b2")]), sd = sqrt(true.2["q"]))
sim.mar2 <- rbind(temp1, temp2)


###################################################
### code chunk number 16: Cs_203_mar2-model
###################################################
Z <- matrix(c(1, 0, 0, 1, 0, 0, 0, 0), 2, 4)
B1 <- matrix(list(0), 2, 2)
diag(B1) <- "b1"
B2 <- matrix(list(0), 2, 2)
diag(B2) <- "b2"
B <- matrix(list(0), 4, 4)
B[1:2, 1:2] <- B1
B[1:2, 3:4] <- B2
B[3:4, 1:2] <- diag(1, 2)
U <- matrix(0, 4, 1)
Q <- matrix(list(0), 4, 4)
Q[1, 1] <- "q"
Q[2, 2] <- "q"
A <- matrix(0, 2, 1)
R <- matrix(0, 2, 2)
pi <- matrix(c(sim.mar2[, 2], sim.mar2[, 1]), 4, 1)
V <- matrix(0, 4, 4)
model.list.2m <- list(Z = Z, B = B, U = U, Q = Q, A = A, 
                      R = R, x0 = pi, V0 = V, tinitx = 1)


###################################################
### code chunk number 19: Cs_204_mar2-fit
###################################################
mar2 <- MARSS(sim.mar2[, 2:TT], model = model.list.2m)


###################################################
### code chunk number 20: Cs_205_mar2-compare
###################################################
model.list.2$x0 <- matrix(sim.mar2[1, 2:1], 2, 1)
mar2a <- MARSS(sim.mar2[1, 2:TT], model = model.list.2)
model.list.2$x0 <- matrix(sim.mar2[2, 2:1], 2, 1)
mar2b <- MARSS(sim.mar2[2, 2:TT], model = model.list.2)


###################################################
### code chunk number 21: Cs_206_compare-mars
###################################################
print(cbind(true = true.2[2:4], est.mar2 = coef(mar2, type = "vector"), est.mar2a = coef(mar2a, type = "vector"), est.mar2b = coef(mar2b, type = "vector")))


###################################################
### code chunk number 23: Cs_301_sim-ar3-data
###################################################
TT <- 100
true.3 <- c(r = 0, b1 = -1.5, b2 = -0.75, b3 = .05, q = 1)
temp3 <- arima.sim(n = TT, list(ar = true.3[c("b1", "b2", "b3")]), 
                   sd = sqrt(true.3["q"]))
sim.ar3 <- matrix(temp3, nrow = 1)


###################################################
### code chunk number 24: Cs_302_set-up-ar3-model
###################################################
Z <- matrix(c(1, 0, 0), 1, 3)
B <- matrix(list("b1", 1, 0, "b2", 0, 1, "b3", 0, 0), 3, 3)
U <- matrix(0, 3, 1)
Q <- matrix(list(0), 3, 3)
Q[1, 1] <- "q"
A <- matrix(0, 1, 1)
R <- matrix(0, 1, 1)
pi <- matrix(sim.ar3[3:1], 3, 1)
V <- matrix(0, 3, 3)
model.list.3 <- list(Z = Z, B = B, U = U, Q = Q, A = A, 
                     R = R, x0 = pi, V0 = V, tinitx = 1)


###################################################
### code chunk number 25: Cs_303_fit-ar3
###################################################
ar3 <- MARSS(sim.ar3[3:TT], model = model.list.3)


###################################################
### code chunk number 26: Cs_304_fit-ar3
###################################################
print(cbind(
  true = true.3[c("b1", "b2", "b3", "q")],
  estimates.no.miss = coef(ar3, type = "vector")
))


###################################################
### code chunk number 27: Cs_401_setseed5
###################################################
set.seed(14)


###################################################
### code chunk number 28: Cs_402_fig-arss-model
###################################################
TT <- 1000 # set long
true.2ss <- c(r = .5, b1 = -1.5, b2 = -0.75, q = .1)
temp <- arima.sim(n = TT, list(ar = true.2ss[c("b1", "b2")]), 
                  sd = sqrt(true.2ss["q"]))
sim.ar <- matrix(temp, nrow = 1)
noise <- rnorm(TT - 1, 0, sqrt(true.2ss["r"]))
noisy.data <- sim.ar[2:TT] + noise


###################################################
### code chunk number 29: Cs_403_fig-arss-model
###################################################
Z <- matrix(c(1, 0), 1, 2)
B <- matrix(list("b1", 1, "b2", 0), 2, 2)
U <- matrix(0, 2, 1)
Q <- matrix(list("q", 0, 0, 0), 2, 2)
A <- matrix(0, 1, 1)
R <- matrix("r")
V <- matrix(0, 2, 2)
pi <- matrix(mean(noisy.data), 2, 1)
model.list.2ss <- list(Z = Z, B = B, U = U, Q = Q, A = A, 
                       R = R, x0 = pi, V0 = V, tinitx = 0)


###################################################
### code chunk number 30: Cs_404_fig-arss-model
###################################################
ar2ss <- MARSS(noisy.data, model = model.list.2ss)


###################################################
### code chunk number 31: Cs_405_fit-arss2-model
###################################################
model.list.2ss.bad <- model.list.2ss
# set R to zero in this model
model.list.2ss.bad$R <- matrix(0)


###################################################
### code chunk number 32: Cs_406_fit-arss2-model
###################################################
ar2ss2 <- MARSS(noisy.data, model = model.list.2ss.bad)


###################################################
### code chunk number 33: Cs_407_fit-arss2-model
###################################################
print(cbind(
  true = true.2ss,
  model.no.error = c(NA, coef(ar2ss2, type = "vector")),
  model.w.error = coef(ar2ss, type = "vector")
))


###################################################
### code chunk number 34: Cs_408_fit-arss2-with-arima
###################################################
arima(noisy.data, order = c(2, 0, 2), include.mean = FALSE)


###################################################
### code chunk number 36: Cs_409_code-to-compare-arss-estimation
###################################################
# This is the code used to make the figure comparing different ways to estimate
# AR parameters from AR with noise data
# sims.exist is just a flag.  Set to FALSE to run the code
if (!exists("sims.exist")) {
  sims.exist <- FALSE
}
if (!sims.exist) {
  nsim <- 200
  TT <- 100
  file <- paste("AR2SS", TT, ".RData", sep = "")
  params <- matrix(0, 8, nsim)
  # sim 2   true.2ss=c(r=.5,b1=0.8,b2=-0.2,q=.1)
  # sim 1
  true.2ss <- c(r = .5, b1 = -1.5, b2 = -0.75, q = .1)

  Z <- matrix(c(1, 0), 1, 2)
  B <- matrix(list("b1", 1, "b2", 0), 2, 2)
  U <- matrix(0, 2, 1)
  Q <- matrix(list("q", 0, 0, 0), 2, 2)
  A <- matrix(0, 1, 1)
  V <- matrix(0, 2, 2)
  R <- matrix("r")
  pi <- matrix(0, 2, 1) # since demeaned
  model.list.2ss <- list(Z = Z, B = B, U = U, Q = Q, A = A, 
                         R = R, x0 = pi, V0 = V, tinitx = 1)

  for (i in 1:nsim) {
    temp <- arima.sim(n = TT, list(ar = true.2ss[2:3]), sd = sqrt(true.2ss[4]))
    sim.ar <- matrix(temp, nrow = 1)
    noise <- rnorm(TT, 0, sqrt(true.2ss[1]))
    noisy.data <- sim.ar + noise
    noisy.data <- as.vector(noisy.data - mean(noisy.data)) # demean
    test.it <- try(arima(noisy.data[2:TT], order = c(2, 0, 2), include.mean = FALSE), silent=TRUE)
    test.it2 <- try(arima(noisy.data[2:TT], order = c(2, 0, 0), include.mean = FALSE), silent=TRUE)     
    while (inherits(test.it, "try-error") | inherits(test.it2, "try-error")) {
      temp <- arima.sim(n = TT, list(ar = true.2ss[2:3]), sd = sqrt(true.2ss[4]))
      sim.ar <- matrix(temp, nrow = 1)
      noise <- rnorm(TT, 0, sqrt(true.2ss[1]))
      noisy.data <- sim.ar + noise
      noisy.data <- as.vector(noisy.data - mean(noisy.data)) # demean
      test.it <- try(arima(noisy.data[2:TT], order = c(2, 0, 2), include.mean = FALSE), silent=TRUE)
      test.it2 <- try(arima(noisy.data[2:TT], order = c(2, 0, 0), include.mean = FALSE), silent=TRUE)
    }
    init.list <- list(Q = matrix(.01, 1, 1), B = matrix(1, 2, 1))
    tmp.kem <- MARSS(noisy.data[2:TT], model = model.list.2ss, inits = init.list, silent = TRUE)
    params[1:2, i] <- coef(tmp.kem)$B
    tmp.bfgs <- MARSS(noisy.data[2:TT], model = model.list.2ss, inits = init.list, silent = TRUE, method = "BFGS")
    # if(any(is.na(tmp.bfgs$states.se)) | any(is.na(tmp.kem$states.se))) print(i)
    params[3:4, i] <- coef(tmp.bfgs)$B
    params[5:6, i] <- test.it$coef[1:2]
    params[7:8, i] <- test.it2$coef[1:2]
    cat(i)
    cat("\n")
    if ((i %% 25) == 0) save(true.2ss, TT, params, file = file)
  }
  sims.exist <- TRUE
  save(true.2ss, TT, params, file = file)
}


###################################################
### code chunk number 37: Cs_410_marssperffig
###################################################
if (sims.exist) { # flag to check that the code to creat the plot has been run
  # This makes a plot of the comparisons
  par(ylbias = -.2, tcl = -.2, cex = .75)
  graphics::boxplot(t(params[c(7, 1, 3, 5, 8, 2, 4, 6), ]), names = c("AR(2)\n", "MARSS\nEM", "MARSS\nBFGS", "ARMA\n(2,2)", "AR(2)\n", "MARSS\nEM", "MARSS\nBFGS", "ARMA\n(2,2)"), ylab = "estimates of the ar coefficients", las = 2)
  points(1:8, apply(params[c(7, 1, 3, 5, 8, 2, 4, 6), ], 1, mean), pch = "x", cex = 1.25)
  par(cex = 1.5)
  axis(side = 3, at = c(2, 6), labels = c(expression(b[1]), expression(b[2])), tick = FALSE, cex = 2)
  lines(c(0, 4.5), c(true.2ss[2], true.2ss[2]), lty = 2)
  lines(c(4.5, 9), c(true.2ss[3], true.2ss[3]), lty = 2)
  abline(v = 4.5)
}


