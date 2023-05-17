###################################################
### code chunk number 26: Cs23_print-B-2
###################################################
# Cleaning up the B matrix for printing
B <- coef(kem.plank.2, type = "matrix")$B[1:4, 1:4]
rownames(B) <- colnames(B) <- c("LP", "SP", "D", "ND")
B[B == 0] <- NA
B.2 <- B
print(B, digits = 2, na.print = "--")


