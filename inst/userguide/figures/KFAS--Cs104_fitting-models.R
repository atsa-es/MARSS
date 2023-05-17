###################################################
### code chunk number 6: Cs104_fitting-models
###################################################
dat <- t(as.matrix(Nile))
rownames(dat) <- "Nile"
fit_em_default <- MARSS(dat, model = mod.nile, silent = TRUE)
inits <- list(Q = matrix(var(Nile)), R = matrix(var(Nile)))
fit_bfgs_default <- MARSS(dat,
  model = mod.nile, inits = inits,
  method = "BFGS", silent = TRUE
)


