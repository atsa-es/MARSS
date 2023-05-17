###################################################
### code chunk number 14: Cs05_fit.model.tinitx1
###################################################
royale.model.4 <- list(
  B = "unconstrained", U = "zero", Q = "diagonal and unequal",
  Z = "identity", R = "zero", tinitx = 1
)
kem.4 <- MARSS(z.royale.dat, model = royale.model.4)


