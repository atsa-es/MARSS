###################################################
### code chunk number 62: Cs704_globaltemp
###################################################
vals <- rbind(
  c(kfas_temp_default$model$Q, kfas_temp_default$model$H[c(1, 2, 4)], -1 * kfas_temp_default$optim.out$value),
  c(coef(marss_temp_default)$Q, coef(marss_temp_default)$R, logLik(marss_temp_default)),
  c(kfas_temp_stoch$model$Q, kfas_temp_stoch$model$H[c(1, 2, 4)], -1 * kfas_temp_stoch$optim.out$value),
  c(coef(marss_temp_stoch_em)$Q, coef(marss_temp_stoch_em)$R, logLik(marss_temp_stoch_em)),
  c(coef(marss_temp_stoch_bfgs)$Q, coef(marss_temp_stoch_bfgs)$R, logLik(marss_temp_stoch_bfgs))
)
rownames(vals) <- c(
  "KFAS default", "MARSS em default",
  "KFAS stoch", "MARSS em stoch", "MARSS bfgs stoch"
)
colnames(vals) <- c("Q", "R1", "Rcov", "R2", "logLik")
round(vals, digits = 5)


