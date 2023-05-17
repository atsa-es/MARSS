###################################################
### code chunk number 6: Cs04_fit.model
###################################################
royale.model.1 <- list(
  Z = "identity", B = "unconstrained",
  Q = "diagonal and unequal", R = "diagonal and unequal",
  U = "zero", tinitx = 1
)
cntl.list <- list(allow.degen = FALSE, maxit = 200)
kem.1 <- MARSS(z.royale.dat, model = royale.model.1, control = cntl.list)


