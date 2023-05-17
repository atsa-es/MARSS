###################################################
### code chunk number 39: Cs27_fit_covar_echo
###################################################
model.list <- list(m = 2, R = "unconstrained")
kemz.temp <- MARSS(dat.spp.1980,
  model = model.list, z.score = TRUE,
  form = "dfa", control = cntl.list, covariates = temp
)
kemz.TP <- MARSS(dat.spp.1980,
  model = model.list, z.score = TRUE,
  form = "dfa", control = cntl.list, covariates = TP
)
kemz.both <- MARSS(dat.spp.1980,
  model = model.list, z.score = TRUE,
  form = "dfa", control = cntl.list, covariates = rbind(temp, TP)
)


