###################################################
### code chunk number 11: Cs_mci_010
###################################################
# create the par list from the output
inits <- coef(kem, what = "par")
bfgs <- MARSS(dat, model = model.list, inits = inits, method = "BFGS")


