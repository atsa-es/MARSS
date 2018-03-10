###############################################################################################################################################
#  glance method for class marssMLE
##############################################################################################################################################
glance.marssMLE <- function (x, ...){
  a = augment(x)
  a = na.omit(a)
  coef.det = cor(a$.fitted,a$y)^2
  ret = data.frame(
    coef.det = cor(a$.fitted,a$y)^2,
    sigma = var(a$.resids),
    df = x$num.params,
    logLik = x$logLik,
    AIC = x$AIC,
    AICc = x$AICc
  )
  if("AICbb"%in%names(x)) ret=cbind(ret, AICbb=x$AICbb)
  if("AICbp"%in%names(x)) ret=cbind(ret, AICbp=x$AICbp)
  if("convergence"%in%names(x)) ret=cbind(ret, convergence=x$convergence) else ret=cbind(ret, convergence=0)
  if("errors"%in%names(x)) ret=cbind(ret, errors=1) else ret=cbind(ret, errors=0)
  ret
}
