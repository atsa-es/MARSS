###################################################
### code chunk number 9: Cs07_fitmodel
###################################################
kem <- MARSS(dat, model = list(
  Z = Z.model,
  Q = Q.model, R = R.model, U = U.model
))


