###################################################
### code chunk number 4: Cs_mci_003
###################################################
U <- x0 <- "zero"
Q <- "unconstrained"
d <- covariates
A <- "zero"
D <- "unconstrained"
R <- "diagonal and equal"
model.list <- list(
  U = U, Q = Q, A = A, R = R,
  D = D, d = d, x0 = x0
)
kem <- MARSS(dat, model = model.list)


