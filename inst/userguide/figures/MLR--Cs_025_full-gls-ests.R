###################################################
### code chunk number 28: Cs_025_full-gls-ests
###################################################
mod4.gls <- gls(Employed ~ 1 + GNP.deflator + GNP + Unemployed
  + Armed.Forces + Population + Year,
correlation = corAR1(), data = longley, method = "ML"
)
mod4.gls.phi <- coef(mod4.gls$modelStruct[[1]], unconstrained = FALSE)
c(mod4.gls.phi, logLik = logLik(mod4.gls))


