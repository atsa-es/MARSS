###################################################
### code chunk number 9: Cs203_structTS-leveltrend
###################################################
mod.list <- list(
  x0 = matrix(c(y[1], 0), 2, 1), U = "zero", tinitx = 0,
  Q = ldiag(c("s2xi", "s2zeta")), R = matrix("s2eps"),
  V0 = matrix(1e+06 * vy, 2, 2) + diag(1e-10, 2), Z = Z, B = B
)
fit3 <- MARSS(as.vector(y), model = mod.list, method = "BFGS")
fit4 <- MARSS(as.vector(y),
  model = mod.list,
  control = list(allow.degen = FALSE)
)


