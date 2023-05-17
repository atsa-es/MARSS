###################################################
### code chunk number 17: Cs2_Code2
###################################################
# Fit the single population model with independent and unequal errors
Z.model <- factor(c(1, 1, 1, 1, 1))
R.model <- "diagonal and unequal"
kem2 <- MARSS(dat, model = list(Z = Z.model, R = R.model))

coef(kem2) # the estimated parameter elements
kem2$logLik # log likelihood
c(kem1$AIC, kem2$AIC) # AICs

plot(kem2, plot.type="model.resids.ytT")


