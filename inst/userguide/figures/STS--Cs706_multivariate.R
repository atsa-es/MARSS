###################################################
### code chunk number 34: Cs706_multivariate
###################################################
Z <- matrix(c(1, 0), n, 2, byrow = TRUE)
mod.list <- c(mod.list.x, mod.list.y, list(Z = Z))
fitm <- MARSS(ym, model = mod.list, method = "BFGS", inits = list(x0 = 0))


