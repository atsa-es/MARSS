###################################################
### code chunk number 25: fitanalysis3
###################################################
model4 <- list(Z = Z4, Q = Q.model, R = R.model4, U = U.model)
kem4 <- MARSS(t(dat4), model = model4)


