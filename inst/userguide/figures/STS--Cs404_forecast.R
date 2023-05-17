###################################################
### code chunk number 19: Cs404_forecast
###################################################
rbind(
  pred1 = fr1$pred, pred2 = fr2$pred$estimate[fr2$ft],
  se1 = fr1$se, se2 = fr2$pred$se[fr2$ft]
)


