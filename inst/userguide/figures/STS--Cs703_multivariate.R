###################################################
### code chunk number 31: Cs703_multivariate
###################################################
vy <- var(y, na.rm = TRUE) / 100
mod.list.x <- list(
  x0 = matrix(list("x0", 0), nrow = 2), tinitx = 1,
  V0 = matrix(1e+06 * vy, 2, 2) + diag(1e-10, 2),
  Q = ldiag(list(q, "qt")),
  B = matrix(c(1, 0, 1, 1), 2, 2),
  U = "zero"
)


