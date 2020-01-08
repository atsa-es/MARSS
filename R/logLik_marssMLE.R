###############################################################################################################################################
#  logLik method for class marssMLE.
##############################################################################################################################################
logLik.marssMLE <- function(object, ...) {
  if (is.null(object$par)) {
    stop("The marssMLE object does not have the par element.  Most likely the model has not been fit.")
  }

  val <- MARSSkf(object, only.logLik = TRUE, return.lag.one = FALSE)$logLik # don't use obj info, recompute ==> object$logLik
  attr(val, "nobs") <- sum(!is.na(object$model$data)) # don't use obj info since user might have updated data ==> object$samp.size
  attr(val, "df") <- length(unlist(object$par)) # on't use obj info since user might have updated model ==> object$num.params
  class(val) <- "logLik"
  val
} # end of logLik.marssMLE
