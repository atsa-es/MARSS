###################################################
### code chunk number 20: Covar_sec6_10_seasonal-fourier-fit
###################################################
model.list <- list(
  B = B, U = U, Q = Q, Z = Z, A = A, R = R,
  C = C, c = c.Four, D = D, d = d
)
seas.mod.3 <- MARSS(dat, model = model.list, control = list(maxit = 1500))


