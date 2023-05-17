###################################################
### code chunk number 7: Cs105_fitting-models
###################################################
mod.nile.stoch <- mod.nile
mod.nile.stoch$x0 <- fit_kfas_stoch$model$a1
mod.nile.stoch$V0 <- fit_kfas_stoch$model$P1
fit_em_stoch <- MARSS(dat, model = mod.nile.stoch, silent = TRUE)
fit_bfgs_stoch <- MARSS(dat,
  model = mod.nile.stoch, inits = inits,
  method = "BFGS", silent = TRUE
)


