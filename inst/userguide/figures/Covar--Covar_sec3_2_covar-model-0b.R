###################################################
### code chunk number 7: Covar_sec3_2_covar-model-0b
###################################################
Q <- "unconstrained"
B <- "diagonal and unequal"
A <- U <- x0 <- "zero"
R <- "diagonal and equal"
d <- covariates
D <- "unconstrained"
y <- dat
model.list <- list(
  B = B, U = U, Q = Q, Z = Z, A = A,
  R = R, D = D, d = d, x0 = x0
)
control.list <- list(maxit = 1500)
kem <- MARSS(y, model = model.list, control = control.list)


