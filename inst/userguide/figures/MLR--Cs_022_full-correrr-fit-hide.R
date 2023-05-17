###################################################
### code chunk number 25: Cs_022_full-correrr-fit-hide
###################################################
inits <- list(A = coef(mod3.em)$A, D = coef(mod3.em)$D)
mod4.em <- MARSS(Employed, model = longley.correrr.model, inits = inits)
mod4.bfgs <- MARSS(Employed,
  model = longley.correrr.model,
  inits = inits, method = "BFGS"
)


