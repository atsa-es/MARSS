###################################################
### code chunk number 34: Cs409_conf-int
###################################################
conf_marss2_t1 <- predict(fit_marss,
  type = "ytt1",
  interval = "confidence", level = 0.95
)
head(conf_marss2_t1$pred)


