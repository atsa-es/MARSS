###################################################
### code chunk number 8: Cs202_structTS-leveltrend
###################################################
vy <- var(y, na.rm = TRUE) / 100
B <- matrix(c(1, 0, 1, 1), 2, 2)
Z <- matrix(c(1, 0), 1, 2)
# fitx parameters at fit1 values
mod.list <- list(
  x0 = matrix(c(y[1], 0), 2, 1), U = "zero", tinitx = 0,
  Q = diag(fit1$coef[1:2]), R = matrix(fit1$coef[3]),
  V0 = matrix(1e+06 * vy, 2, 2), Z = Z, B = B
)
fit2 <- MARSS(as.vector(y),
  model = mod.list, fit = FALSE,
  control = list(trace = -1)
)
fit2$par <- fit2$start # otherwise par is NULL since fit=FALSE


