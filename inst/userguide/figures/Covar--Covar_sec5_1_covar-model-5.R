###################################################
### code chunk number 10: Covar_sec5_1_covar-model-5
###################################################
model.list$R <- diag(0.16, 2)
kem <- MARSS(dat, model = model.list)


