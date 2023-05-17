###################################################
### code chunk number 32: Cs407_conf-int
###################################################
conf_kfas_t1 <- predict(fit_kfas$model,
  interval = "confidence",
  se.fit = TRUE, filtered = TRUE
)
head(conf_kfas_t1)


