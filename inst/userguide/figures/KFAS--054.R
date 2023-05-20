###################################################
### code chunk number 54: KFAS.Rnw:686-693
###################################################
mod.nile.stoch.kfas <- mod.nile.stoch
mod.nile.stoch.kfas$Q <- matrix(fit_kfas_NA$model$Q)
mod.nile.stoch.kfas$R <- matrix(fit_kfas_NA$model$H)
fit_marss_NA <- MARSS(as.vector(NileNA),
  model = mod.nile.stoch.kfas,
  inits = inits, method = "BFGS", silent = TRUE
)


