###################################################
### code chunk number 64: Cs706_globaltemp
###################################################
mod.list <- list(
  Z = matrix(kfas_temp_stoch$model$Z, ncol = 1),
  R = kfas_temp_stoch$model$H[, , 1],
  U = matrix(0),
  A = matrix(0, 2, 1),
  Q = matrix(kfas_temp_stoch$model$Q[, , 1]),
  x0 = kfas_temp_stoch$model$a1,
  V0 = kfas_temp_stoch$model$P1,
  tinitx = 1
)
marss_test <- MARSS(t(GlobalTemp), model = mod.list)


