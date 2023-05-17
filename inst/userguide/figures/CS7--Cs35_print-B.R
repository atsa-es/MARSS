###################################################
### code chunk number 38: Cs35_print-B
###################################################
# Cleaning up the B matrix for printing
B.5 <- coef(kem.plank.5, type = "matrix")$B[1:4, 1:4]
rownames(B.5) <- colnames(B.5) <- c("LP", "SP", "D", "ND")
B.5[B.5 == 0] <- NA
print(B.5, digits = 2, na.print = "--")


