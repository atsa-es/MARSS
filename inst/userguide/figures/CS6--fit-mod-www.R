###################################################
### code chunk number 18: fit-mod-www
###################################################
dat <- as.vector(WWWusage)
kem.3 <- MARSS(dat,
  model = model, inits = list(x0 = matrix(0, 2, 1)),
  control = list(maxit = 20), silent = TRUE
)
kem.3 <- MARSS(dat,
  model = model, inits = kem.3,
  method = "BFGS", silent = TRUE
)
summary(kem.3)


