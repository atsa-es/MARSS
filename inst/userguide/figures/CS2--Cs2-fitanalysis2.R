###################################################
### code chunk number 14: Cs2-fitanalysis2
###################################################
Z.model <- factor(c(1, 1, 1, 1, 1))
R.model <- "diagonal and unequal"
kem2 <- MARSS(dat, model = list(Z = Z.model, R = R.model))


