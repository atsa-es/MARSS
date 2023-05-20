###################################################
### code chunk number 39: fit_covar
###################################################
if (!saved.res) {
  model.list <- list(m = 2, R = "unconstrained")
  kemz.temp <- MARSS(dat.spp.1980,
    model = model.list, z.score = TRUE, form = "dfa",
    control = big.maxit.cntl.list, covariates = temp
  )
  kemz.TP <- MARSS(dat.spp.1980,
    model = model.list, z.score = TRUE, form = "dfa",
    control = big.maxit.cntl.list, covariates = TP
  )
  kemz.both <- MARSS(dat.spp.1980,
    model = model.list, z.score = TRUE, form = "dfa",
    control = big.maxit.cntl.list, covariates = rbind(temp, TP)
  )
} else {
  # reload to re-define kemz.temp etc since prev chunk redefined
  load(paste("./manual_files/", file, sep = ""))
}


