###################################################
### code chunk number 31: Cs406_conf-int
###################################################
pred_marss2 <- predict(fit_marss,
  type = "ytT",
  interval = "prediction", level = 0.95
)


