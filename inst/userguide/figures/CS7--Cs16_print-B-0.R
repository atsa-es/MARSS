###################################################
### code chunk number 19: Cs16_print-B-0
###################################################
# Cleaning up the B matrix for printing
B.0 <- coef(kem.plank.0, type = "matrix")$B[1:4, 1:4]
rownames(B.0) <- colnames(B.0) <- c("LP", "SP", "D", "ND")
print(B.0, digits = 2)


