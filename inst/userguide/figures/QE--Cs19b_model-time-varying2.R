###################################################
### code chunk number 23: Cs19b_model-time-varying2
###################################################
U1 <- matrix(c(rep("t1", 4), "hc"), 5, 1)
U2 <- matrix(c(rep("t2", 4), "hc"), 5, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
Ut[, , 1:floor(TT / 2)] <- U1
kem.tv <- MARSS(dat, model = list(U = Ut, Q = Qde))


