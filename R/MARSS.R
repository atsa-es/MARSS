MARSS <- function(y,
                  model = NULL,
                  inits = NULL,
                  miss.value = as.numeric(NA),
                  method = "kem",
                  form = "marxss",
                  fit = TRUE,
                  silent = FALSE,
                  control = NULL,
                  fun.kf = c("MARSSkfas", "MARSSkfss"),
                  ...) {
  # If user did not pass in fun.kf, then MARSSkfas will be used by default. MARSSkfss will be tried if that fails
  if (missing(fun.kf)) missing.fun.kf <- FALSE else missing.fun.kf <- TRUE
  fun.kf <- match.arg(fun.kf)
  allowed.methods <- get("allowed.methods", envir = pkg_globals)
  # Some error checks depend on an allowable method
  if (!method %in% allowed.methods) {
    stop(paste("method must be one of:", allowed.methods))
  }
  ## Start by checking the data, since if the data have major problems then the rest of the code
  ## will have problems
  if (!missing(miss.value)) {
    stop("miss.value is deprecated in MARSS.  Replace missing values in y with NA.\n")
  }
  if (is.null(y)) {
    stop("MARSS: No data (y) passed in.", call. = FALSE)
  }
  if (!(is.vector(y) | is.matrix(y) | inherits(y, "ts"))) stop("MARSS: Data (y) must be a vector, matrix (time going across columns) or ts/mts object.", call. = FALSE)
  if (length(y) == 0) stop("MARSS: Data (y) is length 0.", call. = FALSE)
  if (is.vector(y)) y <- matrix(y, nrow = 1)
  if (inherits(y, "ts")) {
    model.tsp <- stats::tsp(y)
    y <- t(y)
  } else {
    model.tsp <- c(1, ncol(y), 1)
  }
  attr(y, "model.tsp") <- model.tsp
  if (any(is.nan(y))) cat("MARSS: NaNs in data are being replaced with NAs.  There might be a problem if NaNs shouldn't be in the data.\nNA is the normal missing value designation.\n")
  y[is.na(y)] <- as.numeric(NA)

  if (inherits(model, "marssMLE")) model <- c(coef(model, type = "matrix"), tinitx=model$model$tinitx, diffuse=model$model$diffuse)
  if (inherits(model, "marssMODEL")) model <- c(marssMODEL.to.list(model), tinitx=model$tinitx, diffuse=model$diffuse)

  MARSS.call <- list(data = y, inits = inits, model = model, control = control, method = method, form = form, silent = silent, fit = fit, fun.kf = fun.kf, ...)

  # First make sure specified equation form has a corresponding function to do the conversion to marssMODEL (form=marss) object
  as.marss.fun <- paste("MARSS.", form[1], sep = "")
  tmp <- try(exists(as.marss.fun, mode = "function"), silent = TRUE)
  if (!isTRUE(tmp)) {
    msg <- paste(" MARSS.", form[1], "() function to construct a marssMODEL (form=marss) object does not exist.\n", sep = "")
    cat("\n", "Errors were caught in MARSS \n", msg, sep = "")
    stop("Stopped in MARSS() due to problem(s) with required arguments.\n", call. = FALSE)
  }

  # Build the marssMODEL object from the model argument to MARSS()
  ## The as.marss.fun() call adds the following to MARSS.inputs list:
  ## $marss Translate model strucuture names (shortcuts) into a marssMODEL (form=marss) object put in $marss
  ## $alt.forms with any alternate marssMODEL objects in other forms that might be needed later
  ## a marssMODEL object is a list(data, fixed, free, tinitx, diffuse)
  ## with attributes model.dims, X.names, form, equation
  ## error checking within the function is a good idea though not required
  ## if changes to the control values are wanted these can be set by changing MARSS.inputs$control
  if (silent == 2) cat("Building the marssMODEL object from the model argument to MARSS().\n")
  MARSS.inputs <- eval(call(as.marss.fun, MARSS.call))
  marss.object <- MARSS.inputs$marss
  if (silent == 2) cat(paste("Resulting model has ", attr(marss.object, "model.dims")$x[1], " state processes and ", attr(marss.object, "model.dims")$y[1], " observation processes.\n", sep = ""))

  ## Check that the marssMODEL object output by MARSS.form() is ok
  ## More checking on the control list is done by is.marssMLE() to make sure the MLEobj is ready for fitting
  if (!identical(MARSS.call$control$trace, -1)) { # turn off all error checking
    if (silent == 2) cat(paste("Checking that the marssMODEL object output by MARSS.", form, "() is ok.\n", sep = ""))
    tmp <- is.marssMODEL(marss.object, method = method)
    if (!isTRUE(tmp)) {
      if (!silent || silent == 2) {
        cat(tmp)
        return(marss.object)
      }
      stop("Stopped in MARSS() due to problem(s) with model specification. If silent=FALSE, marssMODEL object (in marss form) will be returned.\n", call. = FALSE)
    } else {
      if (silent == 2) cat("  marssMODEL is ok.\n")
    }
  }

  ## checkMARSSInputs() call does the following:
  ## Check that the user didn't pass in any illegal arguments
  ## and fill in defaults if some params left off
  ## This does not check model since the marssMODEL object is constructed
  ## via the MARSS.form() function above
  if (silent == 2) cat("Running checkMARSSInputs().\n")
  MARSS.inputs <- checkMARSSInputs(MARSS.inputs, silent = FALSE)


  ###########################################################################################################
  ##  MODEL FITTING
  ###########################################################################################################

  ## MLE estimation
  kem.methods <- get("kem.methods", envir = pkg_globals)
  optim.methods <- get("optim.methods", envir = pkg_globals)
  if (method %in% c(kem.methods, optim.methods)) {

    ## Create the marssMLE object

    MLEobj <- list(marss = marss.object, model = MARSS.inputs$model, control = c(MARSS.inputs$control, silent = silent), method = method, fun.kf = fun.kf)
    # Set the call form since that info needed for MARSSinits
    if (MLEobj$control$trace != -1) {
      MLEobj$call <- MARSS.call
    }

    # This is a helper function to set simple inits for a marss MLE model object
    if (silent == 2) cat("Running MARSSinits() to set start conditions.\n")
    MLEobj$start <- MARSSinits(MLEobj, MARSS.inputs$inits)

    class(MLEobj) <- "marssMLE"

    ## Check the marssMLE object
    ## is.marssMLE() calls is.marssMODEL() to check the model,
    ## then checks dimensions of initial value matrices.
    ## it also checks the control list and add defaults if some values are NULL
    if (MLEobj$control$trace != -1) {
      if (silent == 2) cat("Checking the marssMLE object before fitting.\n")
      tmp <- is.marssMLE(MLEobj)
    }

    # if errors, tmp will not be true, it will be error messages
    if (!isTRUE(tmp)) {
      if (!silent || silent == 2) {
        cat(tmp)
        cat(" The incomplete/inconsistent MLE object is being returned.\n")
      }
      cat("Error: Stopped in MARSS() due to marssMLE object incomplete or inconsistent. \nPass in silent=FALSE to see the errors.\n\n")
      MLEobj$convergence <- 2
      return(MLEobj)
    }
    
    # MLEobj is ok.
    
    if (!fit){
      # will be set to 3 if all fixed
      MLEobj$convergence <- -1
    }
    
    # If all parameters fixed. Set convergence=3 whether or not fit=TRUE
    model.is.fixed <- FALSE
    if (all(unlist(lapply(MLEobj[["marss"]][["free"]], is.fixed)))) {
      model.is.fixed <- TRUE
      if (silent == 2) cat("All parameters fixed. No estimation done.\n")
      MLEobj$convergence <- 3
      MLEobj$par <- list()
      for (el in attr(MLEobj[["marss"]], "par.names")) MLEobj[["par"]][[el]] <- matrix(0, 0, 1)
      kf.out <- try(MARSSkf(MLEobj, only.logLik=TRUE), silent = TRUE)
      if (inherits(kf.out, "try-error")){
        MLEobj$convergence <- 53
        MLEobj$logLik <- NA
      }else{
        # rest of the output will be added below
        MLEobj$logLik <- kf.out$logLik
        kf.out <- try(MARSSkf(MLEobj), silent = TRUE)
        if (inherits(kf.out, "try-error")){
          MLEobj$convergence <- 54
        }
      }
    }

    # If fitting is needed, check that Kalman filter/smoother will run
    if (fit && MLEobj$control$trace != -1 && !model.is.fixed) {
      MLEobj.test <- MLEobj
      MLEobj.test$par <- MLEobj$start
      kftest <- try(MARSSkf(MLEobj.test), silent = TRUE)
      if (inherits(kftest, "try-error")) {
        cat("Error: Stopped in MARSS() before fitting because", fun.kf, "stopped.  Something in the model structure prevents the Kalman filter/smoother (KF) from running.\n Try setting fun.kf to use a different KF function (MARSSkfss or MARSSkfas) or use fit=FALSE and check the model you are trying to fit. You can also try trace=1 to get more progress output. You could try trace=-1 to bypass the initial KF check if you are using method='BFGS' and know the logLik function will run.\n")
        MLEobj.test$convergence <- 2
        return(MLEobj.test)
      }
      if (!kftest$ok) {
        cat(kftest$msg)
        cat(paste("Error: Stopped in MARSS() before fitting because", fun.kf, "stopped.  Something in the model structure prevents the Kalman filter or smoother running.\n Try setting fun.kf to use a different KF function (MARSSkfss or MARSSkfas) or use fit=FALSE and check the model you are trying to fit. You can also try trace=1 to get more progress output.\n", kftest$errors, "\n", sep = ""))
        MLEobj.test$convergence <- 2
        return(MLEobj.test)
      }
      MLEobj.test$kf <- kftest
    }
    # Ey is needed for method=kem
    if (fit && !model.is.fixed && MLEobj$control$trace != -1 && MLEobj$method %in% kem.methods) {
      Eytest <- try(MARSShatyt(MLEobj.test), silent = TRUE)
      if (inherits(Eytest, "try-error")) {
        cat("Error: Stopped in MARSS() before fitting because MARSShatyt() stopped.  Something is wrong with the model structure that prevents MARSShatyt() running.\n\n")
        MLEobj.test$convergence <- 2
        MLEobj.test$Ey <- Eytest
        return(MLEobj.test)
      }
      if (!Eytest$ok) {
        cat(Eytest$msg)
        cat("Error: Stopped in MARSS() before fitting because MARSShatyt returned errors.  Something is wrong with the model structure that prevents function running.\n\n")
        MLEobj.test$convergence <- 2
        MLEobj.test$Ey <- Eytest
        return(MLEobj.test)
      }
    }
    
    # fit and not all parameters estimated
    if (fit && !model.is.fixed) {
        if (silent == 2) cat(paste("Fitting model with ", method, ".\n", sep = ""))
        ## Fit and add param estimates to the object
        if (method %in% kem.methods) {
          MLEobj <- try(MARSSkem(MLEobj), silent = TRUE)
          if (inherits(MLEobj, "try-error")) {
            cat(paste("Error: Stopped in MARSS() before fitting because MARSSkem returned errors.  Try control$trace=1 for more information as the reported error may not be helpful. You can also try method='BFGS' if you are seeing a 'chol' error.\n", MLEobj[1], "\n"))
            MLEobj.test$convergence <- 2
            return(MLEobj.test)
          }
        }
        if (method %in% optim.methods) MLEobj <- MARSSoptim(MLEobj)
      }

      ## Add AIC and AICc and coef to the object
      ## Return as long as $par element and there are no errors, but might not be converged
      if ((MLEobj$convergence %in% c(0, 1, 3, 54)) || (MLEobj$convergence %in% c(10, 11) && method %in% kem.methods)) {
        if (silent == 2) cat("Adding AIC and coefficients.\n")
        MLEobj <- MARSSaic(MLEobj)
        MLEobj$coef <- coef(MLEobj, type = "vector")
      }

      ## Add states.se and ytT.se if no errors.  Return kf and Ey if trace>0
      if (MLEobj$convergence %in% c(0, 1, 3) || (MLEobj$convergence %in% c(10, 11) && MLEobj$method %in% kem.methods)) {
        if (silent == 2) cat("Adding states and states.se.\n")
        kf <- MARSSkf(MLEobj) # use function requested by user; default smoother=TRUE
        MLEobj$states <- kf$xtT
        MLEobj$logLik <- kf$logLik
        if (!is.null(kf[["VtT"]])) {
          m <- attr(MLEobj$marss, "model.dims")[["x"]][1]
          TT <- attr(MLEobj$marss, "model.dims")[["data"]][2]
          states.se <- apply(kf[["VtT"]], 3, function(x) takediag(x))
          # KFS 
          if( any(states.se < 0 & states.se > -1*sqrt(.Machine$double.eps)) ){
            states.se[states.se < 0 & states.se > -1*sqrt(.Machine$double.eps)] <- 0
            MLEobj$errors <- c(MLEobj$errors, paste("\nAlert:", MLEobj$fun.kf, "returned negative values close to machine tolerance on diagonal of VtT and the states.se values for these are set to 0. You may want to use", ifelse(MLEobj$fun.kf=='MARSSkfas', 'MARSSkfss()', 'MARSSkfas()'), "to compute the states standard errors (states.se). See MARSSinfo('negVt') for insight.\n"))
          } 
          if( any(states.se < 0) ){
            states.se[states.se < 0] <- NA
            MLEobj$errors <- c(MLEobj$errors, paste("\nAlert:", MLEobj$fun.kf, "returned negative values on diagonal of VtT and the states.se values for these are set to NA. You may want to use", ifelse(MLEobj$fun.kf=='MARSSkfas', 'MARSSkfss()', 'MARSSkfas()'), "to compute the states standard errors (states.se). See MARSSinfo('negVt') for insight.\n"))
          } 
          states.se <- sqrt(states.se)
          if (m == 1) states.se <- matrix(states.se, 1, TT)
          rownames(states.se) <- attr(MLEobj$marss, "X.names")
        } else {
          states.se <- NULL
        }
        MLEobj[["states.se"]] <- states.se
        Ey <- MARSShatyt(MLEobj)
        MLEobj$ytT <- Ey[["ytT"]]
        if (!is.null(Ey[["OtT"]])) {
          n <- attr(MLEobj$marss, "model.dims")[["y"]][1]
          TT <- attr(MLEobj$marss, "model.dims")[["data"]][2]
          if (n == 1) yy <- matrix(Ey[["OtT"]][, , 1:TT], nrow = 1)
          if (n > 1) {
            yy <- apply(Ey[["OtT"]], 3, function(x) {
              takediag(x)
            })
          }
          ytT.se <- sqrt(yy - MLEobj$ytT^2)
          rownames(ytT.se) <- attr(MLEobj$marss, "Y.names")
        } else {
          ytT.se <- NULL
        }
        MLEobj$ytT.se <- ytT.se

        # Return kf and Ey, if trace = 1
        if (MLEobj$control$trace > 0) {
          if (fun.kf == "MARSSkfas") {
            kfss <- try(MARSSkfss(MLEobj, smoother = FALSE), silent = TRUE)
            if (inherits(kfss, "try-error") || !kfss$ok) {
              msg <- "Not available. MARSSkfss() returned error."
              kfss <- list(Innov = msg, Sigma = msg, J = msg, Kt = msg)
            }
          } else {
            kfss <- kf
          }
          MLEobj$kf <- kf # from above will use function requested by user
          MLEobj$Ey <- Ey # from above
          # these are only returned by MARSSkfss
          MLEobj$Innov <- kfss$Innov
          MLEobj$Sigma <- kfss$Sigma
          MLEobj$J <- kfss$J
          MLEobj$Kt <- kfss$Kt
          if (fun.kf == "MARSSkfss") {
            MLEobj$J0 <- kfss$J0
          } else {
            # From kfss smoother so won't be available if fun.kf=MARSSkfas
            # Line above used smoother = FALSE; here smoother = TRUE
            J0 <- try(MARSSkfss(MLEobj), silent = TRUE)
            if (!inherits(J0, "try-error") && J0$ok) MLEobj$J0 <- J0$J0 else MLEobj$J0 <- "Not available. MARSSkfss() smoother returned error."
          }
        }
        # apply X and Y names various X and Y related elements
        MLEobj <- MARSSapplynames(MLEobj)
      }
    # END Adding info to output #################

    if ((!silent || silent == 2) && MLEobj[["convergence"]] %in% c(0, 1, 3, 10, 11, 12, 54)) {
      print(MLEobj)
    }
    if ((!silent || silent == 2) && !(MLEobj[["convergence"]] %in% c(0, 1, 3, 10, 11, 12, 54))) {
      cat(MLEobj$errors)
    } # 3 added since don't print if fit=FALSE
    if ((!silent || silent == 2) && !fit) 
      print(MLEobj$model)

    return(MLEobj)
  } # end MLE methods

  return("method allowed but it's not in kem.methods or optim.methods so marssMLE object was not created")
}
