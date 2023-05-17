###################################################
### code chunk number 33: Cs705_multivariate
###################################################
Z <- matrix(c(1, 0), 1, 2, byrow = TRUE)
mod.list <- c(mod.list.x, mod.list.y, list(Z = Z))
fitu <- MARSS(ym[1, ], model = mod.list, method = "BFGS", inits = list(x0 = 0))


