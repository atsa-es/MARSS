###################################################
### code chunk number 10: Cs08_cov.wolf.moose.model
###################################################
royale.model.3 <- list(
  Z = "identity", B = "unconstrained",
  Q = "diagonal and unequal", R = "zero", U = "zero",
  C = matrix(list(
    0, "Moose win temp", 0, "Moose win precip",
    0, "Moose sum temp"
  ), 2, 3),
  c = z.score.clim.dat
)


