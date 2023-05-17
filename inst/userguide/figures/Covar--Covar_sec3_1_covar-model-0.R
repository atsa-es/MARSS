###################################################
### code chunk number 6: Covar_sec3_1_covar-model-0
###################################################
Q <- U <- x0 <- "zero"
B <- Z <- "identity"
d <- covariates
A <- "zero"
D <- "unconstrained"
y <- dat # to show relationship between dat & the equation
model.list <- list(
  B = B, U = U, Q = Q, Z = Z, A = A,
  D = D, d = d, x0 = x0
)
kem <- MARSS(y, model = model.list)


