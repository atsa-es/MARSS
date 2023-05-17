###################################################
### code chunk number 29: Cs404_conf-int
###################################################
pred_kfas <- predict(fit_kfas$model,
  interval = "prediction", se.fit = TRUE
)
head(pred_kfas)


