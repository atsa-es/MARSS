###################################################
### code chunk number 9: Cs_mci_008
###################################################
inits <- list(Q = matrix(c(1, 0.5, 0.7), 3, 1), D = 1)
kem <- MARSS(dat, model = model.list, inits = inits)


