###################################################
### code chunk number 28: Cs20_kem3_R_unconstrained
###################################################
big.maxit.cntl.list <- list(minit = 200, maxit = 5000, allow.degen = FALSE)
model.list <- list(m = 2, R = "unconstrained")
the.fit <- MARSS(dat.z, model = model.list, form = "dfa", 
                 control = big.maxit.cntl.list)


