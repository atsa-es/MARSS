###################################################
### code chunk number 14: Cs303_structTS-bsm
###################################################
nf <- frequency(y)
vy <- var(y) / 100
B <- makeB(nf)
Z <- matrix(c(1, 0, 1, rep(0, nf - 2L)), 1, nf + 1)

Q <- ldiag(list("s2xi", "s2zeta", "s2w", 0, 0))
R <- matrix("s2eps")
V0 <- matrix(1e+06 * vy, nf + 1, nf + 1) + diag(1e-10, nf + 1)
mod.list <- list(
  x0 = matrix(c(y[1], rep(0, nf)), ncol = 1),
  U = "zero", A = "zero", tinitx = 0,
  Q = Q, R = R, V0 = V0, Z = Z, B = B
)
fit3 <- MARSS(as.vector(y), model = mod.list, method = "BFGS")
fit4 <- MARSS(as.vector(y),
  model = mod.list,
  control = list(allow.degen = FALSE)
)
fit4$kf <- MARSSkfss(fit4)
fit3$kf <- MARSSkfss(fit3)


