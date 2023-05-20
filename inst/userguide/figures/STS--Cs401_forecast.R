###################################################
### code chunk number 16: Cs401_forecast
###################################################
y <- log10(UKgas)
fit1 <- StructTS(y, type = "BSM")

nf <- frequency(y)
vy <- var(y) / 100
B <- makeB(nf) # defined in the BSM section above
Z <- matrix(c(1, 0, 1, rep(0, nf - 2L)), 1, nf + 1)
V0 <- matrix(1e+06 * vy, nf + 1, nf + 1) + diag(1e-10, nf + 1)
mod.list <- list(
  x0 = matrix(c(y[1], rep(0, nf)), ncol = 1), 
  U = "zero", A = "zero", tinitx = 0,
  Q = diag(c(fit1$coef[1:3], 0, 0)), 
  R = matrix(fit1$coef[4]), 
  V0 = V0, Z = Z, B = B
)
fit2 <- MARSS(as.vector(y), model = mod.list)


