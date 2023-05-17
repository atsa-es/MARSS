###################################################
### code chunk number 37: Cs709_multivariate
###################################################
Z <- matrix(c(1, 0), n, 2, byrow = TRUE)
mod.list <- c(mod.list.x, mod.list.y, list(Z = Z, d = covariate))
fitmc <- MARSS(ymc, model = mod.list, method = "BFGS", inits = list(x0 = 0))


