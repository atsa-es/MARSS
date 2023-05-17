###################################################
### code chunk number 16: fit-mod3
###################################################
model <- mod.nile.3
kem.3 <- MARSS(dat,
  model = model, inits = list(x0 = matrix(c(1000, -4), 2, 1)),
  control = list(maxit = 20), silent = TRUE
)
kem.3 <- MARSS(dat,
  model = model, inits = kem.3,
  method = "BFGS", silent = TRUE
)
summary(kem.3)


