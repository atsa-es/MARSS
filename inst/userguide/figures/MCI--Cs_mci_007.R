###################################################
### code chunk number 8: Cs_mci_007
###################################################
inits <- list(Q = matrix(c(1, 0.5, 0.7), 3, 1))
kem <- MARSS(dat, model = model.list, inits = inits)


