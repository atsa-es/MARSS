###################################################
### code chunk number 61: Cs703_globaltemp
###################################################
mod.list <- list(
  Z = matrix(1, 2, 1),
  R = matrix(c("r1", "c", "c", "r2"), 2, 2),
  U = matrix(0),
  A = matrix(0, 2, 1),
  tinitx = 1
)
marss_temp_default <- MARSS(t(GlobalTemp), model = mod.list)
mod.list$x0 <- kfas_temp_stoch$model$a1
mod.list$V0 <- kfas_temp_stoch$model$P1
marss_temp_stoch_em <- MARSS(t(GlobalTemp), model = mod.list)
# use inits from a short run of EM algorithm
inits <- MARSS(t(GlobalTemp),
  model = mod.list, control = list(maxit = 20),
  silent = TRUE
)
marss_temp_stoch_bfgs <- MARSS(t(GlobalTemp),
  model = mod.list,
  inits = inits, method = "BFGS"
)


