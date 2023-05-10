MARSSfit <- function(x, ...) {
  UseMethod("MARSSfit")
}

MARSSfit.default <- function(x, ...){
  stop("MARSS: Something is wrong. Did not find a method for the marssMLE object. Is method set incorrectly?\n There must be a MARSSfit.<method>() function for each method.")
}

MARSSfit.kem <- function(x, ...) MARSSkem(x)

MARSSfit.BFGS <- function(x, ...) MARSSoptim(x)


