# Attaches Hessian, parSigma and parMean to MLEobj
# Computed at the values in MLEobj$par
# For confidence intervals, this should be the MLEs
MARSShessian <- function(MLEobj, method = c("Harvey1989", "fdHess", "optim")) {
  method <- match.arg(method)
  paramvec <- MARSSvectorizeparam(MLEobj)
  if (length(paramvec) == 0) stop("Stopped in MARSShessian(). No estimated parameter elements thus no Hessian.\n", call. = FALSE)

  MLEobj$parMean <- paramvec

  MLEobj$Hessian <- MARSSFisherI(MLEobj, method = method)

  # When inverting the Hessian, need to deal with NAs in the Hessian
  Hess.tmp <- MLEobj$Hessian
  na.diag <- is.na(diag(Hess.tmp))
  if (any(na.diag)) {
    msg <- "MARSShessian: Hessian has NAs due to numerical errors. parSigma returned but some elements will be NA. See MARSSinfo(\"HessianNA\").\n"
    warning(msg)
    MLEobj$errors <- c(MLEobj$errors, msg)
  }
  # set diags with NA to 1 and non-diag to 0 so that Hessian can be inverted
  diag(Hess.tmp)[na.diag] <- 1
  Hess.tmp[is.na(Hess.tmp)] <- 0 #

  hessInv <- try(solve(Hess.tmp), silent = TRUE)
  if (inherits(hessInv, "try-error")) {
    msg <- "MARSShessian: Hessian could not be inverted to compute the parameter var-cov matrix. parSigma set to NULL.  See MARSSinfo(\"HessianNA\").\n"
    warning(msg)
    MLEobj$parSigma <- NULL
    MLEobj$errors <- c(MLEobj$errors, msg)
  } else {
    # Set NAs values to 0
    diag(hessInv)[na.diag] <- NA
    MLEobj$parSigma <- hessInv
  }

  return(MLEobj)
}
