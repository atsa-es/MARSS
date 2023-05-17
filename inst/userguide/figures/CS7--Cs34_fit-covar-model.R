###################################################
### code chunk number 37: Cs34_fit-covar-model
###################################################
plank.model.5 <- plank.model.4
plank.model.5$B <- B
plank.model.5$C <- C
plank.model.5$Q <- Q
plank.model.5$R <- R
kem.plank.5 <- MARSS(d.plank.dat.w.fish, model = plank.model.5)


