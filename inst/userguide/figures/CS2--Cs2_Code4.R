###################################################
### code chunk number 22: Cs2_Code4
###################################################
Z.model <- factor(c(1, 2, 3, 4, 5))
U.model <- "equal"
Q.model <- "diagonal and equal"
R.model <- "diagonal and unequal"
kem <- MARSS(dat, model = list(
  Z = Z.model,
  U = U.model, Q = Q.model, R = R.model
))


