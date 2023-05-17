###################################################
### code chunk number 15: Cs_013_example2-marss-with-inits
###################################################
inits <- list(A = coef(mod1)$A, D = coef(mod1)$D)
mod2 <- MARSS(Employed,
  model = longley.ar1, inits = inits,
  control = list(maxit = 1000)
)
ests.marss <- c(
  b = coef(mod2)$B, alpha = coef(mod2)$A,
  GNP = coef(mod2)$D[1], Population = coef(mod2)$D[2],
  logLik = logLik(mod2)
)


