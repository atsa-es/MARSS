###################################################
### code chunk number 35: Cs410_conf-int
###################################################
pred_kfas_t1 <- predict(fit_kfas$model,
  interval = "prediction",
  se.fit = TRUE, filtered = TRUE
)
head(pred_kfas_t1)


