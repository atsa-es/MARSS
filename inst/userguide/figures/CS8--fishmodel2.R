###################################################
### code chunk number 20: fishmodel2
###################################################
R.model <- "diagonal and equal"
model2 <- list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem2 <- MARSS(fishdat, model = model2)


