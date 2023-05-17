###################################################
### code chunk number 15: Cs13_model-u-NS
###################################################
regions <- list("N", "N", "N", "S", "S")
U <- matrix(regions, 5, 1)
Q <- matrix(list(0), 5, 5)
diag(Q) <- regions
kem <- MARSS(dat, model = list(U = U, Q = Q))


