###################################################
### code chunk number 16: Cs_014_example2-gls
###################################################
library(nlme)
mod2.gls <- gls(Employed ~ GNP + Population,
  correlation = corAR1(), data = longley, method = "ML"
)
mod2.gls.phi <- coef(mod2.gls$modelStruct[[1]], unconstrained = FALSE)
ests.gls <- c(
  b = mod2.gls.phi, alpha = coef(mod2.gls)[1],
  GNP = coef(mod2.gls)[2], Population = coef(mod2.gls)[3],
  logLik = logLik(mod2.gls)
)


