MARSSfit <- function(x, ...) {
  UseMethod("MARSSfit")
}

MARSSfit.default <- function(x, ...) MARSSfit.kem(x, ...)

MARSSfit.kem <- function(x, ...) MARSSkem(x)

MARSSfit.EM_KFAS <- function(x, ...){
  MARSSkem(x)
}

MARSSfit.EM_KFSS <- function(x, ...){
  MARSSkem(x)
}

MARSSfit.BFGS <- function(x, ...){
  MARSSoptim(x)
}

MARSSfit.BFGS_KFAS <- function(x, ...){
  MARSSoptim(x)
}

MARSSfit.BFGS_KFSS <- function(x, ...){
  MARSSoptim(x)
}
