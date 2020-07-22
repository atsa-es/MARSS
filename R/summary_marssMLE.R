###############################################################################################################################################
#  Summary method for class marssMLE.
###############################################################################################################################################

summary.marssMLE <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  model.dims <- attr(object$model, "model.dims")
  n <- model.dims$y[1]
  m <- model.dims$x[1]
  X.names <- attr(object$model, "X.names")
  Y.names <- attr(object$model, "Y.names")
  
  cat(paste("\n",
            "m: ", m, " state process(es) named ", paste(X.names, collapse = " "), "\n",
            "n: ", n, " observation time series named ", paste(Y.names, collapse = " "), "\n\n",
            sep = ""
  ))
  
  if(is.null(object$par) || 
     all(as.logical(lapply(object$model$free, is.fixed)))){
    cat("No estimated parameters.\n")
    }else{
    print(tidy.marssMLE(object, conf.int=FALSE))
    }

  invisible(object)
}
