###################################################
### code chunk number 31: Cs28_print-C
###################################################
# Cleaning up the C matrix for printing
Cmat <- coef(kem.plank.4, type = "matrix")$C[1:4, 1, drop = FALSE]
rownames(Cmat) <- c("LP", "SP", "D", "ND")
Cmat[Cmat == 0] <- NA
print(Cmat, digits = 2, na.print = "--")


