###################################################
### code chunk number 10: fit-data-koop
###################################################
mod.nile.3 <- list(
  Z = matrix(1), A = matrix(0), R = matrix("r"),
  B = matrix(1), U = matrix(0), Q = matrix("q"),
  x0 = matrix("pi"), tinitx = 1, diffuse = TRUE
)
# kem.3.koop=MARSS(dat, model=mod.nile.3,
#  inits=kem.2em$par, method="BFGS")


