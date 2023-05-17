###################################################
### code chunk number 19: fit_data_big_maxit
###################################################
if (!saved.res) {
  big.maxit.cntl.list <- list(minit = 200, maxit = 5000, allow.degen = FALSE)
  kemz.3 <- MARSS(dat.z, model = dfa.model, control = big.maxit.cntl.list)
}


