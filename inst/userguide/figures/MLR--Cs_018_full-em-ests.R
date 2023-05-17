###################################################
### code chunk number 21: Cs_018_full-em-ests
###################################################
par.names <- c("A.intercept", paste("D", eVar.names, sep = "."))
c(coef(mod3.em, type = "vector")[par.names], logLik = mod3.em$logLik)


