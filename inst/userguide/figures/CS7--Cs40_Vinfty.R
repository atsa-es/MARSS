###################################################
### code chunk number 44: Cs40_Vinfty
###################################################
m <- nrow(B)
vecV <- solve(diag(m * m) - kronecker(B, B)) %*% as.vector(Q)
V_inf <- matrix(vecV, nrow = m, ncol = m)


