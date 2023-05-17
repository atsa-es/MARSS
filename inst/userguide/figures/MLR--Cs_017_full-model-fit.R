###################################################
### code chunk number 20: Cs_017_full-model-fit
###################################################
mod3.em <- MARSS(Employed, model = longley.model)
mod3.bfgs <- MARSS(Employed, model = longley.model, method = "BFGS")


