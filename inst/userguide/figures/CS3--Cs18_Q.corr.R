###################################################
### code chunk number 28: Cs18_Q.corr
###################################################
h <- diag(1 / sqrt(diag(Q.unc)))
Q.corr <- h %*% Q.unc %*% h
rownames(Q.corr) <- unique(Z4)
colnames(Q.corr) <- unique(Z4)

Q.corr


