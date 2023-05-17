###################################################
### code chunk number 27: Cs24_fit-plank-model-3
###################################################
# model 3
plank.model.3 <- plank.model.2
plank.model.3$R <- diag(c(.04, .04, .16, .16))
kem.plank.3 <- MARSS(d.plank.dat, model = plank.model.3)


