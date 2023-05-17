###################################################
### code chunk number 55: Cs603_missing-values
###################################################
conf_kfas_NA <-
  predict(fit_kfas_NA$model, interval = "confidence", filtered = FALSE)
conf_marss_NA <-
  predict(fit_marss_NA, interval = "confidence", type = "ytT", level = 0.95)$pred


