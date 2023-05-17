###################################################
### code chunk number 23: Cs20_print-B-1
###################################################
# Cleaning up the B matrix for printing
B <- coef(kem.plank.1, type = "matrix")$B[1:4, 1:4]
rownames(B) <- colnames(B) <- c("LP", "SP", "D", "ND")
B[B == 0] <- NA
B.1 <- B
print(B, digits = 2, na.print = "--")


