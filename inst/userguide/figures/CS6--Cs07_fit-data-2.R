###################################################
### code chunk number 9: Cs07_fit-data-2
###################################################
kem.2em <- MARSS(dat, model = mod.nile.2, silent = TRUE)
kem.2 <- MARSS(dat,
  model = mod.nile.2,
  inits = kem.2em$par, method = "BFGS", silent = TRUE
)
summary(kem.2)


