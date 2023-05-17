###################################################
### code chunk number 26: Cs401_conf-int
###################################################
conf_kfas <- predict(fit_kfas$model,
  interval = "confidence",
  se.fit = TRUE
)
head(conf_kfas)


