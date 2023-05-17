###################################################
### code chunk number 37: Cs412_conf-int
###################################################
pred_marss2_t1 <- predict(fit_marss,
  type = "ytt1",
  interval = "prediction", level = 0.95
)


