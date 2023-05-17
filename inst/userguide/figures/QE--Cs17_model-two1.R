###################################################
### code chunk number 19: Cs17_model-two1
###################################################
Z <- factor(c("N", "N", "N", "S", "S"))
Q <- "diagonal and equal"
U <- "equal"
kem <- MARSS(dat, model = list(Z = Z, Q = Q, U = U))


