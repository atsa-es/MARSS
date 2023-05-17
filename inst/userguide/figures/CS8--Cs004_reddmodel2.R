###################################################
### code chunk number 6: Cs004_reddmodel2
###################################################
model2 <- model1 # model2 is based on model1
model2$R <- "diagonal and unequal"
kem2 <- MARSS(logRedds, model = model2)


