###################################################
### code chunk number 18: Cs16_model-pan2
###################################################
Z <- factor(c(1, 1, 1, 1, 1))
R <- "diagonal and unequal"
kem <- MARSS(dat, model = list(Z = Z, R = R))


