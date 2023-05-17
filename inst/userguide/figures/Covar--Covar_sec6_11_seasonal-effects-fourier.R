###################################################
### code chunk number 21: Covar_sec6_11_seasonal-effects-fourier
###################################################
C.3 <- coef(seas.mod.3, type = "matrix")$C
# The time series of net seasonal effects
seas.3 <- C.3 %*% c.Four[, 1:period]
rownames(seas.3) <- phytos
colnames(seas.3) <- month.abb


