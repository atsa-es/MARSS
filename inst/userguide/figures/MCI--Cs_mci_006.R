###################################################
### code chunk number 7: Cs_mci_006
###################################################
inits <- list(Q = 1)
kem <- MARSS(dat, model = model.list, inits = inits)
# or
inits <- list(Q = matrix(c(1, 0, 1), 3, 1))
kem <- MARSS(dat, model = model.list, inits = inits)


