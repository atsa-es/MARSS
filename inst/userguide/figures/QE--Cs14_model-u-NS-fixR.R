###################################################
### code chunk number 16: Cs14_model-u-NS-fixR
###################################################
regions <- list("N", "N", "N", "S", "S")
U <- matrix(regions, 5, 1)
Q <- matrix(list(0), 5, 5)
diag(Q) <- regions
R <- diag(0.01, 5)
kem <- MARSS(dat, model = list(U = U, Q = Q, R = R))


