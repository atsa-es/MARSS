###################################################
### code chunk number 10: Cs204_structTS-leveltrend
###################################################
fit2$kf <- MARSSkfss(fit2)
fit3$kf <- MARSSkfss(fit3)
fit4$kf <- MARSSkfss(fit4)
data.frame(
  StructTS = fit1$fitted[, 2], fit2 = fit2$kf$xtt[2, ],
  fit.bfgs = fit3$kf$xtt[2, ], fit.em = fit4$kf$xtt[2, ]
)[1:5, ]


