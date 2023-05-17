###################################################
### code chunk number 22: Cs19_fit-plank-model-1
###################################################
plank.model.1 <- plank.model.0
plank.model.1$Q <- "unconstrained"
kem.plank.1 <- MARSS(d.plank.dat, model = plank.model.1)


