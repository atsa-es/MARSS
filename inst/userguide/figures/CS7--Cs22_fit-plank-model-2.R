###################################################
### code chunk number 25: Cs22_fit-plank-model-2
###################################################
# model 2
plank.model.2 <- plank.model.1
plank.model.2$B <- B.2
kem.plank.2 <- MARSS(d.plank.dat, model = plank.model.2)


