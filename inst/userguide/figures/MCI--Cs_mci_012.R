###################################################
### code chunk number 14: Cs_mci_012
###################################################
dat <- t(harborSeal)
dat <- dat[c(2, nrow(dat)), ]
fit1 <- MARSS(dat)
MCinits <- MARSSmcinit(fit1, control = list(numInits = 10))
fit2 <- MARSS(dat, inits = MCinits)


