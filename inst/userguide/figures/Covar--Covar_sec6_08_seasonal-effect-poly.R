###################################################
### code chunk number 18: Covar_sec6_08_seasonal-effect-poly
###################################################
C.2 <- coef(seas.mod.2, type = "matrix")$C
seas.2 <- C.2 %*% month.cov
rownames(seas.2) <- phytos
colnames(seas.2) <- month.abb


