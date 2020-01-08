#######################################################################################################
#   MARSSsimulate function
#   Parametrically simulates from a MARSS parameter list
#   Only works for marss form.  marxss form needs to be converted to marss before this will work.
#######################################################################################################
simulate.marssMLE <- function(MLEobj, tSteps = NULL, nsim = 1, silent = TRUE, miss.loc = NULL) {
  MARSSsimulate(MLEobj, tSteps = tSteps, nsim = nsim, silent = silent, miss.loc = miss.loc)
}
MARSSsimulate <- function(MLEobj, tSteps = NULL, nsim = 1, silent = TRUE, miss.loc = NULL) {
  # tSteps is the number of time steps to do in each bootstrap of the data
  # miss.loc is an optional (n x tSteps x nsim) matrix specifying where to put missing values
  # if miss.loc is the same for all nsim, can pass in dim=c(n, tSteps)
  if (!inherits(MLEobj, "marssMLE") || is.null(MLEobj$par)) {
    stop("Stopped in MARSSsimulate(). The function requires a marssMLE object with the par element.\n", call. = FALSE)
  }
  MARSSobj <- MLEobj[["marss"]]
  n <- dim(MARSSobj$fixed$A)[1]
  m <- dim(MARSSobj$fixed$x0)[1]

  ###### Error-checking on MARSSobj
  tmp <- is.marssMODEL(MARSSobj, method = MLEobj[["method"]])
  if (!isTRUE(tmp)) {
    if (!silent) cat(tmp)
    stop("Stopped in MARSSsimulate() due to problem with MARSSobj.\n", call. = FALSE)
  }

  ###### Error-checking on the arguments
  msg <- NULL
  if (!is.numeric(nsim)) {
    msg <- c(msg, "Non-numeric nsim argument(s).\n")
  }
  if (!is.null(tSteps) & !is.numeric(tSteps)) {
    msg <- c(msg, "Non-numeric tSteps argument(s).\n")
  }
  if (length(tSteps) > 1 || length(nsim) > 1) {
    msg <- c(msg, "tSteps and nsim must be length 1.\n")
  }
  if (is.numeric(tSteps) && (tSteps != trunc(tSteps) || tSteps <= 0)) {
    msg <- c(msg, "tSteps must be a positive non-zero integer.\n")
  }
  if (is.numeric(nsim) && (nsim != trunc(nsim) || nsim <= 0)) {
    msg <- c(msg, "nsim must be a positive non-zero integer.\n")
  }
  if (!is.null(miss.loc) && (!isTRUE(all.equal(dim(miss.loc), c(n, tSteps)))
  && !isTRUE(all.equal(dim(miss.loc), c(n, tSteps, nsim))))) {
    msg <- c(msg, "Incorrect input arg: miss.loc dim must be n x tSteps or n x tSteps x nsim.\n")
  }

  # Check that if any fixed or free are time-varying, TT = tSteps
  en <- names(MARSSobj$fixed)
  Tmax <- 1
  for (elem in en) {
    Tmax <- max(Tmax, dim(MARSSobj$fixed[[elem]])[3])
    Tmax <- max(Tmax, dim(MARSSobj$free[[elem]])[3])
  }
  if (is.null(tSteps)) tSteps <- Tmax
  if (!(Tmax == 1 || Tmax == tSteps)) {
    msg <- c(msg, "If any fixed or free matrices are time-varying, the dim of time must equal tSteps.\n")
  }

  if (!is.null(msg)) {
    cat("\nErrors were caught in MARSSsimulate\n", msg)
    stop("Stopped in MARSSsimulate() due to argument problem(s).\n", call. = FALSE)
  }
  if (!is.null(miss.loc)) {
    miss.loc.TF <- is.na(miss.loc)
  }

  ##### Set holders for output
  # if user passed in miss.loc dim=c(n,tSteps) assume they wanted that repeated for all nsim's
  # set up miss.loc if user didn't pass it in;
  if (is.null(miss.loc)) { # means no missing values
    miss.loc.TF <- array(FALSE, dim = c(n, tSteps))
  } # 1 means not missing
  if (length(dim(miss.loc.TF)) == 2) miss.loc.TF <- array(miss.loc.TF, dim = c(n, tSteps, nsim))
  # sim.data = array(NA,dim=c(tSteps,n,nsim))
  # sim.states = array(NA,dim=c(tSteps,m,nsim))
  sim.data <- array(as.numeric(NA), dim = c(n, tSteps, nsim))
  sim.states <- array(as.numeric(NA), dim = c(m, tSteps, nsim))

  ##### Set up the progress bar
  drawProgressBar <- FALSE # If the time library is not installed, no prog bar
  if (!silent) { # then we can draw a progress bar
    prev <- progressBar()
    drawProgressBar <- TRUE
  }

  ##### Set up holders
  newData <- matrix(NA, n, tSteps + 1)
  newStates <- matrix(NA, m, tSteps + 1) # States = years x subpops

  #### make a list of time-varying parameters
  time.varying <- list()
  for (elem in names(MARSSobj$free)) {
    if ((dim(MARSSobj$free[[elem]])[3] == 1) & (dim(MARSSobj$fixed[[elem]])[3] == 1)) {
      time.varying[[elem]] <- FALSE
    } else {
      time.varying[[elem]] <- TRUE
    } # not time-varying
  }
  # set the parameters at t=1
  par1 <- parmat(MLEobj, t = 1)

  ##### Construct needed permutation matrices when there are 0s on diag of var-cov matrices
  ##### It is required that the 0 locations are time invariant in Q, R and V0; checked in is.marssMLE()
  Omg1 <- t.Omg1 <- n.not0 <- Omg0 <- t.Omg0 <- list()
  for (elem in c("Q", "R", "V0")) {
    dim.par <- sqrt(dim(MARSSobj$fixed[[elem]])[1])
    Omg1[[elem]] <- t.Omg1[[elem]] <- Omg0[[elem]] <- t.Omg0[[elem]] <- array(0, dim = c(dim.par, dim.par))
    n.not0[[elem]] <- c()
    the.par <- par1[[elem]]
    diag.par <- diag(the.par)
    n.not0[[elem]] <- sum(diag.par != 0)
    I.mat <- diag(1, dim.par)
    if (n.not0[[elem]] == dim.par) {
      Omg1[[elem]] <- t.Omg1[[elem]] <- I.mat
      Omg0[[elem]] <- t.Omg0[[elem]] <- matrix(0, dim.par, dim.par)
    } else {
      if (n.not0[[elem]] == 0) {
        Omg0[[elem]] <- t.Omg0[[elem]] <- I.mat
        Omg1[[elem]] <- t.Omg1[[elem]] <- matrix(0, dim.par, dim.par)
      } else {
        Omg1[[elem]] <- I.mat[diag.par != 0, , drop = FALSE]
        Omg0[[elem]] <- I.mat[diag.par == 0, , drop = FALSE]
        t.Omg1[[elem]] <- t(Omg1[[elem]])
        t.Omg0[[elem]] <- t(Omg0[[elem]])
      }
    }
  } # for over elem

  x0.mat <- par1$x0
  pari <- par1
  for (i in 1:nsim) {
    newStates[, 1] <- x0.mat
    if (n.not0$V0 != 0) { # by def length 1
      V0.mat <- Omg1$V0 %*% par1$V0 %*% t.Omg1$V0
      # rmvnorm returns a 1 x m matrix even if mean is m x 1
      x0.new <- array(rmvnorm(1, mean = Omg1$V0 %*% x0.mat, sigma = V0.mat, method = "chol"), dim = dim(x0.mat))
      newStates[, 1] <- t.Omg1$V0 %*% x0.new + t.Omg0$V0 %*% Omg0$V0 %*% x0.mat
    } else {
      newStates[, 1] <- x0.mat
    }
    for (j in 2:(tSteps + 1)) { # j=1 is t=0, j=2 is t=1
      if (time.varying$R) pari$R <- parmat(MLEobj, "R", t = (j - 1))$R
      # create matrices for observation error
      if (n.not0$R[min(j - 1, length(n.not0$R))] != 0) { # some nonzeros; minus 1 since j=2 is t=1
        R.mat <- Omg1$R %*% pari$R %*% t.Omg1$R
        # rmvnorm returns a T x p matrix and we need p x T
        obs.error <- t(rmvnorm(1, mean = rep(0, n.not0$R), sigma = R.mat, method = "chol"))
        obs.error <- t.Omg1$R %*% obs.error
      } else {
        obs.error <- matrix(0, n, 1)
      } # R all zero
      # create a matrices for process error
      if (n.not0$Q != 0) {
        if (time.varying$Q) pari$Q <- parmat(MLEobj, "Q", t = (j - 1))$Q
        Q.mat <- Omg1$Q %*% pari$Q %*% t.Omg1$Q
        # rmvnorm returns a 1 x p matrix and we need p x 1
        pro.error <- t(rmvnorm(1, mean = rep(0, n.not0$Q), sigma = Q.mat, method = "chol"))
        pro.error <- t.Omg1$Q %*% pro.error
      } else {
        pro.error <- matrix(0, m, 1)
      } # Q all zero
      if (j == 2 && MARSSobj$tinitx == 1) {
        newStates[, 2] <- newStates[, 1]
      } else {
        if (time.varying$B) pari$B <- parmat(MLEobj, "B", t = (j - 1))$B
        if (time.varying$U) pari$U <- parmat(MLEobj, "U", t = (j - 1))$U
        if (time.varying$Z) pari$Z <- parmat(MLEobj, "Z", t = (j - 1))$Z
        if (time.varying$A) pari$A <- parmat(MLEobj, "A", t = (j - 1))$A
        newStates[, j] <- pari$B %*% newStates[, j - 1] + pari$U + pro.error
      }
      newData[, j] <- pari$Z %*% newStates[, j] + pari$A + obs.error
      # odd indexing since the parmat assume you are passin in actual t, but the newStates indexing is 1 off
      # to see it works, sub transform j back to t; say j=3, X(j) is X(2) and X(j-1) is X(1), so we have
      # X(2)=B(2)X(1)+U(2)+pro.error(2)
    }
    newData[miss.loc.TF[, , i]] <- as.numeric(NA)
    newStates <- newStates[, 2:(tSteps + 1)] # make indexing t=1:TT again
    newData <- newData[, 2:(tSteps + 1)]
    sim.data[, , i] <- as.matrix(newData)
    sim.states[, , i] <- as.matrix(newStates)
    # reset newStates and newData to their original dims
    newStates <- matrix(NA, m, tSteps + 1)
    newData <- matrix(NA, n, tSteps + 1)
    # Draw the progress bar if silent=F and time library is installed
    if (drawProgressBar) {
      prev <- progressBar(i / nsim, prev)
    }
  } # end of for loop for nsim

  return(list(sim.states = sim.states, sim.data = sim.data, MLEobj = MLEobj, miss.loc = miss.loc, tSteps = tSteps, nsim = nsim))
}
