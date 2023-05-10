#######################################################################################################
#   MARSSboot
#   This creates bootstrap parameter estimates
#   This is an MLE function and uses marssMLE objects
#   return(list(boot.params=boot.params, boot.data=boot.data, model=MLEobj, nboot=nboot, output=output, sim=sim, param.gen=param.gen, control=control))
#######################################################################################################
MARSSboot <- function(MLEobj, nboot = 1000, output = "parameters", sim = "parametric",
                      param.gen = "MLE", control = NULL, silent = FALSE) {

  # all equation numbers refer to Chapter 6, Shumway & Stoffer
  # tSteps is the number of time steps to do in each bootstrap of the data
  #     the default is to take this from the number of years in the original data
  # nboot is the number of bootstrap replicates to do
  # MLEobj is a marssMLE object (needs the structure and parameter elements of the list)
  # output tells what output to produce; this can be a vector
  #      data : returns the simulated (or bootstrapped) data.
  #      parameters : return the ML parameter estimates from the bootstrapped data
  #      all : return both
  # sim = the type of bootstrapping to do to create boot.data (if output included "data").
  #      Choices are
  #      parametric uses parametric bootstrapping
  #      innovations uses Stoffer and Walls method
  # param.gen = how to generate the parameters.
  #      MLE --- generate sim data and estimate params via MLE
  #      hessian  -- get params from estimated hessian matrix
  # control is a list which holds options for the estimation function (see help file)
  # silent controls whether a progress bar is shown

  ###### Error-checking on the arguments
  msg <- NULL
  if (!is.numeric(nboot)) {
    msg <- c(msg, " Incorrect function arg. nboot must be numeric.\n")
  }
  if (is.numeric(nboot) && any(nboot != trunc(nboot), nboot < 0, length(nboot) != 1)) {
    msg <- c(msg, " Incorrect function arg. nboot must be a positive non-zero integer.\n")
  }
  if (FALSE %in% (output %in% c("all", "data", "parameters"))) {
    msg <- c(msg, " Incorrect function arg. output must be either all, data, parameters (in quotes).\n")
  }
  if (FALSE %in% (sim %in% c("parametric", "innovations"))) {
    msg <- c(msg, " Incorrect function arg. sim must be either parametric or innovations (see help file).\n")
  }
  if (FALSE %in% (param.gen %in% c("MLE", "hessian"))) {
    msg <- c(msg, " Incorrect function arg. est.method must be either MLE or hessian (see help file).\n")
  }
  if (FALSE %in% (silent %in% c(TRUE, FALSE))) {
    msg <- c(msg, " Incorrect function arg. silent must be TRUE or FALSE")
  }
  if (!(MLEobj$method %in% allowed.methods)) {
    msg <- c(msg, paste0(" MLE object method must be in ", paste0(allowed.methods, collapse = ", "), "."))
  }
  if (!is.null(msg)) {
    cat("\n", "Errors were caught in MARSSboot \n", msg, sep = "")
    stop("Stopped in MARSSboot() due to problem(s) with function arguments.\n", call. = FALSE)
  }

  # mle.object$control is control args for estimation routine
  # control here is just if you want to change the control values in mle.object$control
  # code replaces mle.object$control element with the one in the passed in control obj
  if (!is.null(control)) {
    if (!is.list(control)) {
      msg <- " Incorrect MARSSbootstrap arg. control must be a list.\n"
      cat("\n", "Errors were caught in MARSSboot \n", msg, sep = "")
      stop("Stopped in MARSSboot() due to problem(s) with function arguments.\n", call. = FALSE)
    }
    EM.null <- NULL
    en <- c("minit", "maxit", "abstol", "min.degen.iter", "degen.lim", "min.iter.conv.test", "conv.test.deltaT", "conv.test.slope.tol")
    for (el in en) {
      null.flag <- (is.null(control[[el]]))

      if (!null.flag) {
        null.flag <- (!is.numeric(control[[el]]) || control[[el]] <= 0)
        if (!null.flag) MLEobj$control[[el]] <- control[[el]] # if control$el passed in and it is ok, reset MLEobj$control
      }
      EM.null <- c(EM.null, null.flag)
    }
    problem <- any(EM.null)
    if (problem) {
      msg <- paste(" Invalid control", en[EM.null], "\n")
      cat("\n", "Errors were caught in MARSSboot \n", msg, sep = "")
      stop("Stopped in MARSSboot() due to problem(s) with control argument.\n", call. = FALSE)
    }
    ##  trace must be set to (0) FALSE otherwise too much memory would be consumed
    MLEobj$control$trace <- 0
    MLEobj$control$silent <- TRUE # prevent output from the estimation function
  }

  # Check for marssMLE properness; check that MLEobj$marss$data exists
  if (!inherits(MLEobj, "marssMLE")) {
    stop("Stopped in MARSSboot(). Object must be class marssMLE.\n", call. = FALSE)
  }

  # Check that it has par added on
  if (is.null(MLEobj[["par"]])) {
    msg <- " MLE object is missing parameter estimates.\n"
    cat("\n", "Errors were caught in MARSSboot \n", msg, sep = "")
    stop("Stopped in MARSSboot() due to problem(s) with MLE object.\n", call. = FALSE)
  }

  if (output == "all") output <- c("data", "parameters")

  # Check that if param.gen=hessian, the user hasn't set output to just "data"
  if (all(param.gen == "hessian", "data" %in% output)) {
    if (!MLEobj$control$silent) warning("MARSSboot: when param.gen=hessian, no boot.data are output")
    output <- output[output != "data"]
    if (length(output) == 0) {
      msg <- "  When param.gen=hessian, output arg must be set to parameters.\n"
      cat("\n", "Errors were caught in MARSSboot \n", msg, sep = "")
      stop("Stopped in MARSSboot() due to problem(s) with function arguments.\n", call. = FALSE)
    }
  }
  if (param.gen == "hessian") sim <- "none"

  # A little renamimg for code readability
  marss.model <- MLEobj[["marss"]]
  ##### Now check for any inconsistency in the passed in arguments
  if (any(is.na(marss.model$data)) & sim == "innovations") {
    msg <- "  Innovations bootstrapping uses the innovations resampling and can only be done if there are no missing values in the data.\n"
    cat("\n", "Errors were caught in MARSSboot \n", msg, sep = "")
    stop("Stopped in MARSSboot() due to problem(s) with function arguments.\n", call. = FALSE)
  }
  ##### Set holders for output
  boot.params <- NA
  boot.data <- NA # dummy values
  # these are the arrays for output
  if ("data" %in% output) {
    boot.data <- array(NA, dim = c(dim(as.matrix(marss.model$data)), nboot))
    rownames(boot.data) <- attr(marss.model[["data"]], "Y.names")
  }
  if ("parameters" %in% output) {
    boot.params <- array(NA, dim = c(length(MARSSvectorizeparam(MLEobj)), nboot))
    rownames(boot.params) <- names(MARSSvectorizeparam(MLEobj))
  }

  ##### Read in model parameters
  marss.model <- MLEobj[["marss"]]
  par.dims <- attr(marss.model, "model.dims")
  m <- par.dims[["x"]][1]
  TT <- par.dims[["data"]][2] # length of time series
  n <- par.dims[["y"]][1]

  ##### If using hessian to generate boot parameters, check if parSigma is already set
  if (param.gen == "hessian" & (is.null(MLEobj[["parSigma"]]) | is.null(MLEobj[["parMean"]]))) {
    if (!silent) cat("MARSSboot: Computing the Hessian.  This might take awhile.\n")
    MLEobj.hessian <- MARSShessian(MLEobj) # returns a transformed MLEobj with hessian; parSigma; and parMean matrix onto the mle.model object
    if (is.null(MLEobj.hessian[["parSigma"]])) { # parSigma is not set if hessian cannot be inverted
      stop("Stopped in MARSSboot() because Hessian could not be inverted to estimate the parameter var-cov matrix", call. = FALSE)
    }
  }

  ##### Set up the progress bar
  drawProgressBar <- FALSE # If the time library is not installed, no prog bar
  if (!silent) { # then we can draw a progress bar
    prev <- progressBar()
    drawProgressBar <- TRUE
  }

  ##### This part creates and stores the bootstrap data (not parameters yet just data)
  if (sim == "parametric") {
    boot.data <- MARSSsimulate(MLEobj, miss.loc = marss.model$data, tSteps = TT, nsim = nboot)$sim.data
  } # make new data

  if (sim == "innovations") {
    boot.data <- MARSSinnovationsboot(MLEobj, nboot = nboot)$boot.data
  }

  ##### This part generates the bootstrap parameter
  if ("parameters" %in% output) {
    boot.control <- MLEobj$control
    boot.control$silent <- TRUE # turn off output from MARSSkem()
    mle.object <- MLEobj
    mle.object$control <- boot.control
    for (i in 1:nboot) {
      if (param.gen == "MLE") {
        newmod <- marss.model # marss form
        newmod[["data"]] <- array(boot.data[, , i], dim = dim(boot.data)[1:2], dimnames = dimnames(marss.model$data))
        mle.object[["marss"]] <- newmod # we are resetting the marss object
        # fitting function determined by the MLEobj$method
        boot.model <- MARSSfit(mle.object)
        boot.params[, i] <- MARSSvectorizeparam(boot.model)
      } # if MLE

      if (param.gen == "hessian") {
        if (any(is.na(MLEobj.hessian$parSigma))) stop("Stopped in MARSSboot(): The variance-covariance matrix (parSigma) returned by MARSShessian() has NAs.  See MARSSinfo(\"HessianNA\").\n", call. = FALSE)
        hess.params <- rmvnorm(1, mean = MLEobj.hessian$parMean, sigma = MLEobj.hessian$parSigma, method = "chol")
        boot.params[, i] <- hess.params
      } # if hessian

      # Draw the progress bar if silent=FALSE and time library is installed
      if (drawProgressBar) prev <- progressBar(i / nboot, prev)
    } # end nboot loop
  } # end if parameters in output

  return(list(boot.params = boot.params, boot.data = boot.data, marss = marss.model, nboot = nboot, output = output, sim = sim, param.gen = param.gen))
}
