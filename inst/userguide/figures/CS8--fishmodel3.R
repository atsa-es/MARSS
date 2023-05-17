###################################################
### code chunk number 21: fishmodel3
###################################################
R.model <- matrix(list(0), 9, 9)
diag(R.model) <- list("1", "2", "4", "4", "4", "4", "4", "5", "6")
model3 <- list(Z = Z.model, Q = Q.model, R = R.model, U = U.model)
kem3 <- MARSS(fishdat, model = model3)


