#' Generic for fitting MARSS models
#' 
#' Uses the method of a marssMLE class object. Will call a function such as [MARSSkem()], [MARSSoptim()] or `MARSStmb()` in the {marssTMB} package.
#' @param x a [marssMLE] object.
#' @param ... additional arguments for the fitting function
MARSSfit <- function(x, ...) {
  UseMethod("MARSSfit")
}

MARSSfit.default <- function(x, ...){
  stop("MARSS: Something is wrong. Did not find a method for the marssMLE object. Is method set incorrectly?\n There must be a MARSSfit.<method>() function for each method.")
}

MARSSfit.kem <- function(x, ...) MARSSkem(x)

MARSSfit.BFGS <- function(x, ...) MARSSoptim(x)


