###################################################
### code chunk number 9: Cs107_fitting-models
###################################################
vals <- rbind(
  c(fit_kfas_default$model$Q, fit_kfas_default$model$H, -1 * fit_kfas_default$optim.out$value),
  c(coef(fit_em_default)$Q, coef(fit_em_default)$R, logLik(fit_em_default)),
  c(coef(fit_bfgs_default)$Q, coef(fit_bfgs_default)$R, logLik(fit_bfgs_default)),
  c(fit_kfas_stoch$model$Q, fit_kfas_stoch$model$H, -1 * fit_kfas_stoch$optim.out$value),
  c(coef(fit_em_stoch)$Q, coef(fit_em_stoch)$R, logLik(fit_em_stoch)),
  c(coef(fit_bfgs_stoch)$Q, coef(fit_bfgs_stoch)$R, logLik(fit_bfgs_stoch)),
  c(fit_marss_kfas$model$Q, fit_marss_kfas$model$H, -1 * fit_marss_kfas$optim.out$value)
)
rownames(vals) <- c(
  "KFAS default", "MARSS em default", "MARSS bfgs default",
  "KFAS stoch", "MARSS em stoch", "MARSS bfgs stoch", "KFAS w marss kfas model"
)
colnames(vals) <- c("Q", "R", "logLik")
vals


