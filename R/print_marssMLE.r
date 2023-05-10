###############################################################################################################################################
#  Print method for class marssMLE.
##############################################################################################################################################
print.marssMLE <- function(x, digits = max(3, getOption("digits") - 4), ..., what = "fit", form = NULL, silent = FALSE) {
  # load needed package globals
  kem.methods <- get("kem.methods", envir = pkg_globals)
  optim.methods <- get("optim.methods", envir = pkg_globals)
  nlminb.methods <- get("nlminb.methods", envir = pkg_globals)
  
  # Check that x has a marssMODEL object
  if (!inherits(x$model, "marssMODEL")) {
    cat("\nThe model element of your marssMLE object is not class marssMODEL.\n")
    cat("Are you using a marssMLE object from a MARSS version before 3.5?\n")
    cat("Type MARSSinfo(\"modelclass\") for help.\n")
    return()
  }
  # Check that x has a marssMODEL object
  if (!inherits(x$marss, "marssMODEL")) {
    cat("\nThe marss element of your marssMLE object is not class marssMODEL.\n")
    cat("Are you using a marssMLE object from a MARSS version before 3.5?\n")
    cat("Type MARSSinfo(\"modelclass\") for help.\n")
    return()
  }

  #### Error-checking that 'what' can be printed

  # No output for these if model was not fit
  what.vals <- c("fit", "start", "inits", "par", "logLik", "paramvector", "par.se", "par.bias", "par.upCI", "par.lowCI", "xtT", "states", "ytT", "states.se", "model.residuals", "state.residuals", "kfs", "Ey", "states.cis", names(x$model$fixed))
  if (is.null(x$par) && any(what %in% what.vals)) {
    if (silent) {
      return()
    }
    cat("marssMLE object $par element is NULL.  Parameters have not been estimated.\n")
    if (x$convergence == 2) {
      cat("ERROR: The fitting function returned a try-error (crashed).\n marssMLE object $convergence arg = 2\n")
    }
    if (x$convergence == 52 || x$convergence == 62) {
      cat(
        "WARNING: Estimation was stopped due to errors.\n",
        "see $errors element of marssMLE object to view errors\n"
      )
    }
    if (x$convergence == 53) {
      cat(
        "ERROR: Estimation was stopped in optim() due to errors returned by MARSSkf.\n",
        "No parameter estimates are available.\n marssMLE object $convergence arg = 53\n"
      )
    }
    return()
  }

  # No params estimated
  what.vals <- c("paramvector", "start", "inits", "par.se", "par.bias", "par.upCI", "par.lowCI")
  if (any(what %in% what.vals) && all(unlist(lapply(x$model$free, is.fixed)))) {
    if (!silent) cat("No estimated parameters so no paramvector, parameters CIs or standard errors.\n")
    return()
  }

  # Set up error message
  conv54msg <- "No parameter CIs with Hessian, Kalman filter/smoother output or residuals calculations possible because MARSSkf (the Kalman filter/smoother) returns an error with the fitted model. Try MARSSinfo('optimerror54') for insight.\n"
  what.vals <- c("par.se", "par.bias", "par.upCI", "par.lowCI", "xtT", "states", "ytT", "states.se", "model.residuals", "state.residuals", "kfs", "Ey", "states.cis")
  if (any(what %in% what.vals) && identical(x$convergence, 54)) {
    if (!silent) cat(conv54msg)
    return()
  }
  ##### End error-checking on 'what'

  # Set up what to return
  return.obj <- list()
  what.to.print <- what
  orig.x <- x

  ## Run the print_form function to put x$par in the correct form
  if (is.null(form)) { # allow user to print using a different form than x$model form
    form <- attr(x[["model"]], "form")
  }
  print.fun <- paste("print_", form[1], sep = "")
  tmp <- try(exists(print.fun, mode = "function"), silent = TRUE)
  if (isTRUE(tmp)) {
    # the print function can print or return an updated x to use for printing
    # it changes the $par element to match $model instead of $marss
    x <- eval(call(print.fun, x))
  } else {
    x <- print_null(x)
  } # if print_form is missing use print_marss

  if (identical(what.to.print, "fit")) {
    if (silent) {
      return(orig.x)
    } else {
      invisible(orig.x)
    }

    ## all parameters were fixed; don't use convergence since this won't catch fixed and conv=54
    if (all(unlist(lapply(x$model$free, is.fixed)))) {
      cat("\nAll parameters were fixed so none estimated.\n")
      cat("Log-likelihood:", x$logLik, "\n")
      if (!is.null(x[["AIC"]])) cat("AIC:", x$AIC, "  ")
      if (!is.null(x[["AICc"]])) cat("AICc:", x$AICc, "  ")
      if (!is.null(x[["AICbb"]])) cat("AICbb(innov):", x$AICbb, "  ")
      if (!is.null(x[["AICbp"]])) cat("AICbp(param):", x$AICbp, "  ")
      cat("\n \n")
      if (identical(x$convergence, 54)) {
        cat("WARNING: MARSSkf (the Kalman filter/smoother) returned an error with the model. No states estimates.\n")
      }

      invisible(orig.x)
    } else {
      cat("\nMARSS fit is\n")
      cat(paste("Estimation method:", x$method, "\n"))
      if (x$method %in% kem.methods) {
        cat(paste("Convergence test: conv.test.slope.tol = ", x$control$conv.test.slope.tol, ", abstol = ", x$control$abstol, "\n", sep = ""))
      }
      if (x$convergence == 0) {
        if (x$method %in% c(optim.methods, nlminb.methods)) {
          cat("Estimation converged in", x$numIter, "iterations. \n")
        } else {
          if (x$numIter > x$control$minit) {
            cat("Estimation converged in", x$numIter, "iterations. \n")
          } else {
            cat("Algorithm ran", x$numIter, "(=minit) iterations and convergence was reached. \n")
          }
        }
      }
      if (x$convergence == 1) {
        if (x$method %in% kem.methods) tmp.msg <- paste("Neither abstol nor log-log convergence test were passed.\n", sep = "") else tmp.msg <- x$iter.record$message
        if(x$method %in% c(kem.methods, optim.methods) || (x$method %in% nlminb.methods & x$numIter >= x$control$maxit)) cat("WARNING: maxit reached at ", x$control$maxit, " iter before convergence.\n", tmp.msg,
          "The likelihood and params are not at the MLE values.\n",
          "Try setting control$maxit higher.\n"
        )
        if(x$method %in% nlminb.methods && x$numIter < x$control$maxit)
          cat("WARNING: nlminb() ran for", x$numIter, "iterations and stopped before maxit with a convergence warning.\nTreat the parameter values and logLik with caution. Message:", tmp.msg, ".\n")
      }
      if (x$convergence == 2) cat("Invalid MLE object. Try-error returned by fitting function. \n")
      if (x$convergence == 3) cat("All parameters fixed. No convergence information. \n")
      if (x$method %in% kem.methods && x$convergence == 10) {
        tmp.msg <- paste("maxit (=", x$control$maxit, ") reached before log-log convergence.\n", sep = "")
        cat(
          "WARNING: Abstol convergence only no log-log convergence.\n", tmp.msg,
          "The likelihood and params might not be at the ML values.\n",
          "Try setting control$maxit higher.\n"
        )
      }
      if (x$method %in% kem.methods && x$convergence == 12) {
        tmp.msg <- paste("No log-log convergence info because maxit (=", x$control$maxit, ") < min.iter.conv.test (=", x$control$min.iter.conv.test, ").\n", sep = "")
        cat(
          "WARNING: Abstol convergence only no info on log-log convergnece.\n", tmp.msg,
          "The likelihood and params might not be at the ML values.\n",
          "Try setting control$maxit higher.\n"
        )
      }
      if (x$method %in% kem.methods && x$convergence == 11) {
        tmp.msg <- paste("maxit (=", x$control$maxit, ") reached before abstol convergence.\n", sep = "")
        cat(
          "WARNING: log-log convergence only no abstol convergence.\n", tmp.msg,
          "The likelihood and params might not be at the ML values.\n",
          "Try setting control$maxit higher.\n"
        )
      }
      if (x$method %in% optim.methods && x$convergence == 10) {
        cat("WARNING: degeneracy of the Nelder-Mead simplex. \n")
      }
      if (x$method %in% kem.methods && x$convergence %in% c(52, 62)) {
        cat(
          "WARNING: Estimation was stopped at iteration", x$numIter, "due to errors.\n",
          "Parameter estimates are those at the last iteration before stopping\n",
          "see $errors element of marssMLE object to view errors\n"
        )
      }
      if (x$method %in% optim.methods && (x$convergence == 51 || x$convergence == 52)) {
        cat("WARNING: error or warning from the L-BFGS-B method; see component message for details\n")
      }
      if (!(x$convergence %in% c(0, 1, 2, 3, 10, 11, 12, 13, 51, 52, 62))) {
        cat("WARNING: convergence test errors\n")
      }
      if ((x$convergence %in% c(54))) {
        cat("WARNING: Kalman filter-smoother error.\n")
      }
      cat("Log-likelihood:", x$logLik, "\n")
      if (!is.null(x[["AIC"]])) cat("AIC:", x$AIC, "  ")
      if (!is.null(x[["AICc"]])) cat("AICc:", x$AICc, "  ")
      if (!is.null(x[["AICbb"]])) cat("AICbb(innov):", x$AICbb, "  ")
      if (!is.null(x[["AICbp"]])) cat("AICbp(param):", x$AICbp, "  ")
      cat("\n \n")

      # set up the paramvector
      paramvector <- coef(orig.x, type = "vector")

      cmat <- as.matrix(paramvector)
      colnames(cmat) <- "Estimate"
      if (!is.null(x[["par.se"]])) {
        cmat <- cbind(cmat, coef(orig.x, type = "vector", what = "par.se"))
        colnames(cmat) <- c("ML.Est", "Std.Err")
      }

      if (!is.null(x[["par.lowCI"]]) & !is.null(x[["par.upCI"]])) {
        cmat <- cbind(cmat, coef(orig.x, type = "vector", what = "par.lowCI"))
        colnames(cmat)[dim(cmat)[2]] <- "low.CI"
        cmat <- cbind(cmat, coef(orig.x, type = "vector", what = "par.upCI"))
        colnames(cmat)[dim(cmat)[2]] <- "up.CI"
      }
      if (!is.null(x[["par.bias"]])) {
        par.bias <- coef(orig.x, type = "vector", what = "par.bias")
        cmat <- cbind(cmat, par.bias)
        colnames(cmat)[dim(cmat)[2]] <- "Est.Bias"
        cmat <- cbind(cmat, as.matrix(paramvector) + par.bias)
        colnames(cmat)[dim(cmat)[2]] <- "Unbias.Est"
      }
      # printCoefmat(cmat, digits=digits)
      print(cmat, digits = digits)
      if (x$model$tinitx == 1) {
        cat("Initial states (x0) defined at t=1\n")
      } else {
        cat("Initial states (x0) defined at t=0\n")
      }

      cat("\n")
      if (is.null(x[["par.se"]])) cat("Standard errors have not been calculated. \n")
      if (!is.null(x[["par.lowCI"]]) & !is.null(x[["par.upCI"]])) {
        cat(paste("CIs calculated at alpha = ", x$par.CI.info$alpha, " via method=", x$par.CI.info$method, sep = ""), "\n")
        if (x$par.CI.info$method == "hessian" && (any(is.na(x[["par.lowCI"]])) || any(is.na(x[["par.upCI"]])) | any(is.na(x[["par.se"]])))) {
          cat("There are NAs in the Hessian matrix. Type MARSSinfo(\"HessianNA\") for details.\n")
        }
      } else {
        cat("Use MARSSparamCIs to compute CIs and bias estimates.\n")
      }
      if (!is.null(x[["par.bias"]])) {
        cat(paste("Bias calculated via", x$par.CI.method, "bootstrapping with", x$par.CI.nboot, "bootstraps. \n"))
      }

      if (!is.null(x[["errors"]])) {
        if (length(x$errors) > 10) {
          cat(x$errors[1])
          cat(paste(length(x$errors) - 1, " ", ifelse(grepl("warnings", x$errors[1]), "warnings", "errors"), ". First 10 shown.  Type cat(object$errors) to see the full list.\n ", sep = ""))
          cat(x$errors[2:11])
        } else {
          cat(x$errors)
        }
      }
    }
    cat("\n")
    invisible(orig.x)
  } else { # what!=fit

    for (what in what.to.print) {
      if (what == "fit") return.obj[[what]] <- orig.x
      if (what == "model") {
        # note at beginning of this function x$model was replaced with the model
        # in call$form form
        if (!silent) print.marssMODEL(x$model)
        return.obj[[what]] <- x$model
      }
      if (what == "par") {
        if (!silent) {
          cat("List of the estimated values in each parameter matrix\n")
          print(x$par)
          cat("\n")
        }
        return.obj[[what]] <- x$par
      }
      if (what == "logLik") {
        if (!silent) {
          cat(paste("Log-likelihood: ", x$logLik, sep = ""))
          cat("\n")
        }
        return.obj[[what]] <- x$logLik
      }
      if (what %in% c("par.se", "par.bias", "par.upCI", "par.lowCI")) {
        if (is.null(x[[what]])) {
          if (what != "par.bias") {
            cat(paste("No ", what, " element on your marssMLE object.\nRunning MARSSparamCIs with method=hessian\n", sep = ""))
            orig.x <- MARSSparamCIs(orig.x)
          } else {
            cat("No par.bias element on your marssMLE object.\nRunning MARSSparamCIs with method=parametric to get parameter bias estimate.\n")
            orig.x <- MARSSparamCIs(x, method = "parametric")
          }
        }
        if (!silent) {
          cat(paste("Estimated ", what, ":\n"))
          print(coef(orig.x, type = "vector", what = what))
          cat("\n")
        }
        return.obj[[what]] <- orig.x[[what]]
      }
      if (what == "paramvector") {
        paramvector <- coef(orig.x, type = "vector")
        if (!silent) {
          cat("vector of the estimated parameter values\n")
          print(paramvector)
          cat("\n")
        }
        return.obj[[what]] <- paramvector
      }
      if (what == "xtT" | what == "states") {
        if (!silent) {
          cat("Smoothed state estimates\n")
          print(x$states)
          cat("\n")
        }
        return.obj[[what]] <- orig.x$states
      }
      if (what == "data") {
        if (!silent) {
          model.tsp <- attr(orig.x$model$data, "model.tsp")
          attr(orig.x$model$data, "model.tsp") <- NULL
          print(ts(t(orig.x$model$data), start = model.tsp[1], frequency = model.tsp[3]))
        }
        return.obj[[what]] <- orig.x$model$data
      }
      if (what == "ytT") {
        if (is.null(x[["ytT"]])) {
          ytT <- MARSShatyt(orig.x)$ytT
        } else {
          ytT <- x$ytT
        }
        if (!silent) {
          cat("Expected value of Y(t) conditioned on all the data.  If data are present ytT=y.\n")
          print(ytT)
          cat("\n")
        }
        return.obj[[what]] <- ytT
      }
      if (what == "states.se") {
        if (!silent) {
          cat("Smoothed state standard errors.\n")
          print(orig.x$states.se)
          cat("\n")
        }
        return.obj[[what]] <- orig.x$states.se
      }
      if (what == "model.residuals") {
        res <- MARSSresiduals(orig.x, type = "tt1")$model.residuals
        if (!silent) {
          cat("Innovations model residuals.\n")
          print(res)
          cat("\n")
        }
        return.obj[[what]] <- res
      }
      if (what == "state.residuals") {
        res <- MARSSresiduals(orig.x, type = "tT")$state.residuals
        if (!silent) {
          cat("Smoothed state residuals.\n")
          print(res)
          cat("\n")
        }
        return.obj[[what]] <- res
      }
      if (what == "kfs") {
        kf <- MARSSkf(orig.x)
        if (x$fun.kf == "MARSSkfas") {
          tmp <- try(MARSSkfss(orig.x), silent = TRUE) # MARSSkfas doesn't return Innov or Sigma
          if (!inherits(tmp, "try-error")) {
            # will be NULL if kfss failed
            kf$Innov <- tmp$Innov
            kf$Sigma <- tmp$Sigma
          }
        }
        if (!silent) {
          cat("Kalman filter and smoother output is verbose.  Assign print output and look at assigned object.\n")
          cat("\n")
        }
        return.obj[[what]] <- kf
      }
      if (what == "Ey") {
        Ey <- MARSShatyt(orig.x)
        if (!silent) {
          cat("Expectations involving y conditioned on all the data.  See ?MARSShatyt. \n Output is verbose.  Assign print output and look at assigned object.\n")
          cat("\n")
        }
        return.obj[[what]] <- Ey
      }
      if (what == "states.cis") {
        if (!silent) cat("Approximate 95% confidence intervals for the states using 1.95*states.se.\n\n")
        rpt.list <- list(up95CI = orig.x$states + qnorm(0.975) * orig.x$states.se, est = orig.x$states, low95CI = orig.x$states - qnorm(0.975) * orig.x$states.se)
        if (!silent) {
          print(rpt.list)
          cat("\n")
        }
        return.obj[[what]] <- rpt.list
      }
      if (what == "start" | what == "inits") {
        if (!silent) {
          cat("Initial parameter values used for optimization.\n\n")
          print(coef(orig.x, type = "vector", what = "start"))
          cat("\n")
        }
        return.obj[[what]] <- coef(orig.x, type = "vector", what = "start")
      }
      fixed <- x$model$fixed
      free <- x$model$free
      # if there is info on the parameter dimensions in the model list, use it
      par.dims <- attr(x$model, "model.dims")
      if (what %in% names(fixed)) {
        par.dim <- par.dims[[what]]
        the.par <- parmat(x, what, t = 1:max(dim(fixed[[what]])[3], dim(free[[what]])[3]), dims = par.dim, model.loc = "model")[[what]]
        if (!silent) {
          cat(paste("Parameter matrix ", what, "\n", sep = ""))
          print(the.par)
          cat("\n")
        }
        return.obj[[what]] <- the.par
      }
      if (what %in% names(x$model) && !(what == "data")) {
        if (!silent) {
          cat(paste("Model list element ", what, "\n", sep = ""))
          print(x$model[[what]])
          cat("\n")
        }
        return.obj[[what]] <- x$model[[what]]
      }
    } # for what in what.list
    if (length(return.obj) == 0) return.obj <- NULL
    if (length(return.obj) == 1) return.obj <- return.obj[[1]]
    if (silent) {
      return(return.obj)
    } else {
      invisible(return.obj)
    }
  } # what is not identical to "fit"
} # end of print.marssMLE

print_null <- function(x) {
  return(x)
}
