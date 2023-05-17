###################################################
### code chunk number 8: Cs06_print-wolf.B
###################################################
wolf.B <- coef(kem.2, type = "matrix")$B
rownames(wolf.B) <- colnames(wolf.B) <- rownames(royale.dat)
print(wolf.B, digits = 2)


