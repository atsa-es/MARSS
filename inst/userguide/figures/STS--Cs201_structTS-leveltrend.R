###################################################
### code chunk number 7: Cs201_structTS-leveltrend
###################################################
y <- log10(forecast:::subset.ts(UKgas, quarter = 2))
fit1 <- StructTS(y, type = "trend")


