#######################################################################################################
#   MARSShessian.numerical functions
#
#   Numerical estimation of the Hessian of the negative log-likelihood function
#     at the parameters values in MLEobj$par.  For the Hessian to be used
#     for confidence interval calculation, this should be the MLEs.
#     Two functions are available for the Hessian calculation, fdHess() and optim().
#
#   Adds Hessian, parameter var-cov matrix, and parameter mean to a marssMLE object
#######################################################################################################
MARSShessian.numerical <- function(MLEobj, fun = "fdHess") {
  ##
  kfNLL <- function(x, MLEobj = NULL) { # NULL assignment needed for optim call syntax
    MLEobj <- MARSSvectorizeparam(MLEobj, x)
    negLL <- MARSSkf(MLEobj, only.logLik = TRUE, return.lag.one = FALSE)$logLik
    -1 * negLL
  }

  paramvector <- MARSSvectorizeparam(MLEobj)

  # Hessian and gradient
  if (fun == "fdHess") {
    out <- fdHess(paramvector, function(paramvector, MLEobj) kfNLL(paramvector, MLEobj), MLEobj)
    Hessian <- out$Hessian
  }
  if (fun == "optim") {
    # maxit set to 0 so that the Hessian is computed at the values in the MLEobj
    out <- optim(paramvector, kfNLL, MLEobj = MLEobj, method = "BFGS", hessian = TRUE, control = list(maxit = 0))
    Hessian <- out$hessian
  }

  rownames(Hessian) <- names(paramvector)
  colnames(Hessian) <- names(paramvector)
  # MLEobj$gradient = out$gradient

  return(Hessian)
}
