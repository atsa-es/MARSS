###################################################
### code chunk number 8: Cs106_fitting-models
###################################################
marss_kfas_model <- MARSSkfas(fit_em_stoch,
  return.kfas.model = TRUE,
  return.lag.one = FALSE
)$kfas.model
marss_kfas_model$Q[1, 1, 1] <- NA
marss_kfas_model$H[1, 1, 1] <- NA
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_marss_kfas <- fitSSM(marss_kfas_model, kinits, method = "BFGS")


