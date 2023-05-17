###################################################
### code chunk number 5: Cs103_structTS-level
###################################################
fit2$kf <- MARSSkfss(fit2)
fit3$kf <- MARSSkfss(fit3)
fit4$kf <- MARSSkfss(fit4)
df <- data.frame(
  StructTS = fit1$fitted, fit2 = fit2$kf$xtt[1, ],
  fit.bfgs = fit3$kf$xtt[1, ], fit.em = fit4$kf$xtt[1, ]
)
head(df)


