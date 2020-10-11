MARSSresiduals <- function(object, ..., type = c("tT", "tt1", "tt"), normalize = FALSE, silent = FALSE, fun.kf = c("MARSSkfas", "MARSSkfss")) {
  type <- match.arg(type)
  fun.kf <- match.arg(fun.kf)
  if (is.null(object[["par"]])) {
    stop("MARSSresiduals: The marssMLE object does not have the par element.  Most likely the model has not been fit.", call. = FALSE)
  }
  if (object[["convergence"]] == 54) {
    stop("MARSSresiduals: optim() successfully fit this model but MARSSkf (the Kalman filter/smoother) returns an error with the fitted model. Try MARSSinfo('optimerror54') for insight.", call. = FALSE)
  }
  
  if (type == "tT") val <- MARSSresiduals.tT(object, ..., normalize = normalize, silent = silent, fun.kf = fun.kf)
  if (type == "tt1") val <- MARSSresiduals.tt1(object, ..., normalize = normalize, silent = silent, fun.kf = fun.kf)
  if (type == "tt") val <- MARSSresiduals.tt(object, ..., normalize = normalize, silent = silent, fun.kf = fun.kf)
  return(val)
}
