###################################################
### code chunk number 28: Cs403_conf-int
###################################################
conf_marss2 <- predict(fit_marss,
  type = "ytT",
  interval = "confidence", level = 0.95
)
head(conf_marss2$pred)


