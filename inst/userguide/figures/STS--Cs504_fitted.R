###################################################
### code chunk number 24: Cs504_fitted
###################################################
fitted3 <- MARSSkfss(fit2)$xtt
fitted3 <- ts(t(fitted3[1:3, ]))
plot(fitted3)


