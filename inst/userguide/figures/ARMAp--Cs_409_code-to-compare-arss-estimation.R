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
  model.list.2ss <- list(
    Z = Z, B = B, U = U, Q = Q, A = A,
    R = R, x0 = pi, V0 = V, tinitx = 1
  )

  for (i in 1:nsim) {
    temp <- arima.sim(n = TT, list(ar = true.2ss[2:3]), sd = sqrt(true.2ss[4]))
    sim.ar <- matrix(temp, nrow = 1)
    noise <- rnorm(TT, 0, sqrt(true.2ss[1]))
    noisy.data <- sim.ar + noise
    noisy.data <- as.vector(noisy.data - mean(noisy.data)) # demean
    test.it <- try(arima(noisy.data[2:TT], order = c(2, 0, 2), include.mean = FALSE), silent = TRUE)
    test.it2 <- try(arima(noisy.data[2:TT], order = c(2, 0, 0), include.mean = FALSE), silent = TRUE)
    while (inherits(test.it, "try-error") | inherits(test.it2, "try-error")) {
      temp <- arima.sim(n = TT, list(ar = true.2ss[2:3]), sd = sqrt(true.2ss[4]))
      sim.ar <- matrix(temp, nrow = 1)
      noise <- rnorm(TT, 0, sqrt(true.2ss[1]))
      noisy.data <- sim.ar + noise
      noisy.data <- as.vector(noisy.data - mean(noisy.data)) # demean
      test.it <- try(arima(noisy.data[2:TT], order = c(2, 0, 2), include.mean = FALSE), silent = TRUE)
      test.it2 <- try(arima(noisy.data[2:TT], order = c(2, 0, 0), include.mean = FALSE), silent = TRUE)
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


