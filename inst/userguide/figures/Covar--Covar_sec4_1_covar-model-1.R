###################################################
### code chunk number 8: Covar_sec4_1_covar-model-1
###################################################
R <- A <- U <- "zero"
B <- Z <- "identity"
Q <- "equalvarcov"
C <- "unconstrained"
x <- dat # to show the relation between dat & the equations
model.list <- list(
  B = B, U = U, Q = Q, Z = Z, A = A,
  R = R, C = C, c = covariates
)
kem <- MARSS(x, model = model.list)


