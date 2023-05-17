###################################################
### code chunk number 13: Cs_mci_011
###################################################
######################################################################################################   MARSSmcinit function
#   Does a simple MonteCarlo initialization of the EM routine
#   The function uses a number of MARSS utility functions accessed with MARSS:::
#######################################################################################################
MARSSmcinit <- function(MLEobj,
                        control = list(
                          numInits = 500, numInitSteps = 10,
                          MCbounds = list(
                            B = c(0, 1), U = c(-1, 1), Q = c(1, 1),
                            Z = c(0, 1), A = c(-1, 1), R = c(1, 1), x0 = c(-1, 1)
                          )
                        ),
                        silent = FALSE) {
  control.default <- list(numInits = 500, numInitSteps = 10, MCbounds = list(B = c(0, 1), U = c(-1, 1), Q = c(1, 1), Z = c(0, 1), A = c(-1, 1), R = c(1, 1), x0 = c(-1, 1)))
  if (!is.null(control)) {
    if (!is.list(control)) stop("MARSSmcinit: control must be a list")
    if (any(!(names(control) %in% names(control.default)))) stop(paste("MARSSmcinit: allowed control list elements are", names(control.default)))
    control.new <- control.default
    for (i in names(control)) control.new[[i]] <- control[[i]]
    control <- control.new
  }
  drawProgressBar <- FALSE
  if (!silent) { # then we can draw a progress bar
    cat("\n")
    cat("> Starting Monte Carlo Initializations\n")
    prev <- MARSS:::progressBar() # this is an internal function to MARSS
    drawProgressBar <- TRUE
  }
  MODELobj <- MLEobj[["marss"]]
  y <- MODELobj$data
  par.dims <- attr(MODELobj, "model.dims")
  m <- par.dims[["x"]][1]
  n <- par.dims[["y"]][1]
  TT <- par.dims[["data"]][2]
  ## YM matrix for handling missing values
  YM <- matrix(as.numeric(!is.na(y)), n, TT)
  # Make sure the missing vals in y are zeroed out
  y[YM == 0] <- 0

  free.tmp <- MODELobj$free
  dim.tmp <- list(Z = c(n, m), A = c(n, 1), R = c(n, n), B = c(m, m), U = c(m, 1), Q = c(m, m), x0 = c(m, 1))
  bounds.tmp <- control$MCbounds
  init <- bestinits <- MLEobj$start
  bestLL <- -1.0e10

  # loop over numInits: # of random draws of initial values
  for (loop in 1:control$numInits) {
    init.loop <- init

    # Draw random values
    en <- c("Z", "A", "R", "B", "U", "Q", "x0")
    for (el in en) {
      dim.param <- dim.tmp[[el]]
      if (!MARSS:::is.fixed(free.tmp[[el]])) { # is.fixed is a utility func in MARSS
        bounds.param <- bounds.tmp[[el]]
        # use the first fixed and free in a temporally varying model; arbitrary
        tmp <- list(f = MARSS:::sub3D(MODELobj$fixed[[el]], t = 1), D = MARSS:::sub3D(MODELobj$free[[el]], t = 1))
        if (el %in% c("Q", "R")) { # random starts drawn from a wishart dist
          if (bounds.param[1] < dim.param[1]) {
            df <- dim.param[1]
          } else {
            df <- bounds.param[1]
          }
          S <- diag(bounds.param[2], dim.param[1])
          # draw a random matrix from wishart
          tmp.random <- MARSS:::rwishart(df, S) / df
          # reapply the sharing and fixed constraints
          par.random <- solve(t(tmp$D) %*% tmp$D) %*% t(tmp$D) %*% (MARSS:::vec(tmp.random) - tmp$f)
        } else {
          par.random <- matrix(runif(dim(tmp$D)[2], bounds.param[1], bounds.param[2]), dim(tmp$D)[2], 1)
          if (el %in% c("B")) {
            tmp.max <- max(abs(eigen(par.random, only.values = TRUE)$values))
            # rescale to bring the max abs eigenvalues to between 0 and 1
            par.random <- par.random / (tmp.max / runif(1, .01, .99))
          }
          if (el %in% c("x0")) {
            x0init <- init$x0 # where the original start is
            x.lo <- ifelse(x0init > 0, exp(bounds.param[1]) * x0init, exp(bounds.param[2]) * x0init)
            x.hi <- ifelse(x0init > 0, exp(bounds.param[2]) * x0init, exp(bounds.param[1]) * x0init)
            par.random <- matrix(runif(dim(tmp$D)[2], x.lo, x.hi), dim(tmp$D)[2], 1)
          }
        }
      } else {
        par.random <- matrix(0, 0, 1)
      }
      init.loop[[el]] <- par.random
    }

    ## Call MARSSkem() with these inits
    MLEobj$start <- init.loop
    MLEobj$control$maxit <- control$numInitSteps
    MLEobj$control$minit <- 1
    MLEobj$control$silent <- TRUE # don't output
    MLEobj <- MARSSkem(MLEobj) # get new fit using this init

    if (drawProgressBar == TRUE) prev <- MARSS:::progressBar(loop / control$numInits, prev)

    ## Check whether the likelihood is the best observed
    ## Only use bootstrap param draws where loglike did not go down during numInitSteps
    if (MLEobj$logLik > bestLL) {
      # update the best initial parameter estimates
      bestinits <- MLEobj$par
      bestLL <- MLEobj$logLik
    }
  } # end numInits loop

  return(bestinits)
}


