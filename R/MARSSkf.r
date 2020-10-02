#######################################################################################################
#   MARSSkf function
#   Utility function to choose the Kalman filter and smoother
#######################################################################################################
MARSSkf <- function(MLEobj, only.logLik = FALSE, return.lag.one = TRUE, return.kfas.model = FALSE, newdata = NULL, smoother = TRUE) {
  if (is.null(MLEobj$par)) {
    stop("Stopped in MARSSkf(): par element of marssMLE object is required.\n")
  }
  if (!missing(newdata)) {
    data.dim <- attr(MLEobj[["model"]], "model.dims")[["data"]]
    if (!identical(dim(newdata), data.dim)) {
      stop("Stopped in MARSSkf(): newdata must be the same dimensions as the data used to fit the model.\n")
    }
    if (!is.numeric(newdata)) {
      stop("Stopped in MARSSkf(): newdata must be numeric.\n")
    }
    if (dim(MLEobj[["model"]][["free"]][["x0"]])[2] != 0) {
      message("MARSSkf(): x0 was estimated and this x0 will be used with newdata.\n")
    }
    MLEobj[["model"]][["data"]] <- newdata
    MLEobj[["marss"]][["data"]] <- newdata
  }
  if (MLEobj$fun.kf == "MARSSkfss") {
    return(MARSSkfss(MLEobj, smoother = smoother))
  }
  if (MLEobj$fun.kf == "MARSSkfas") {
    return(MARSSkfas(MLEobj, only.logLik = only.logLik, return.lag.one = return.lag.one, return.kfas.model = return.kfas.model))
  }
  return(list(ok = FALSE, errors = "kf.function does not specify a valid Kalman filter and smoother function."))
}
