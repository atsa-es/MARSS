###################################################
### code chunk number 10: Cs_mci_009
###################################################
inits <- list(D = coef(kem, what = "par")$D)
kem <- MARSS(dat, model = model.list, inits = inits)


