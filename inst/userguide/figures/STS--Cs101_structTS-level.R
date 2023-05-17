###################################################
### code chunk number 3: Cs101_structTS-level
###################################################
y <- window(treering, start = 0, end = 20)

fit1 <- StructTS(y, type = "level")


