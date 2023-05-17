###################################################
### code chunk number 53: Cs602_missing-values
###################################################
rbind(
  MARSS = c(
    Q = coef(fit_marss_NA, type = "matrix")$Q,
    R = coef(fit_marss_NA, type = "matrix")$R,
    logLik = logLik(fit_marss_NA)
  ),
  KFAS = c(
    Q = fit_kfas_NA$model$Q,
    R = fit_kfas_NA$model$H,
    logLik = -1 * fit_kfas_NA$optim.out$value
  )
)


