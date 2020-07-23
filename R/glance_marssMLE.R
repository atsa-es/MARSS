###############################################################################################################################################
#  glance method for class marssMLE. For package broom
##############################################################################################################################################
glance.marssMLE <- function(x, ...) {
  a <- residuals.marssMLE(x)
  a <- subset(a, a$type=="model")
  a <- na.omit(a)
  ret <- data.frame(
    coef.det = cor(a$.pred, a$value)^2,
    sigma = var(a$.resids),
    df = x$num.params,
    logLik = x$logLik,
    AIC = x$AIC,
    AICc = x$AICc,
    stringsAsFactors = FALSE
  )
  if ("AICbb" %in% names(x)) ret <- cbind(ret, AICbb = x$AICbb)
  if ("AICbp" %in% names(x)) ret <- cbind(ret, AICbp = x$AICbp)
  if ("convergence" %in% names(x)) ret <- cbind(ret, convergence = x$convergence) else ret <- cbind(ret, convergence = 0)
  if ("errors" %in% names(x)) ret <- cbind(ret, errors = 1) else ret <- cbind(ret, errors = 0)
  ret
}
