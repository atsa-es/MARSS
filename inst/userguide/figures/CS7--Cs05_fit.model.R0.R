###################################################
### code chunk number 7: Cs05_fit.model.R0
###################################################
royale.model.2 <- list(
  Z = "identity", B = "unconstrained",
  Q = "diagonal and unequal", R = "zero", U = "zero"
)
kem.2 <- MARSS(z.royale.dat, model = royale.model.2)


