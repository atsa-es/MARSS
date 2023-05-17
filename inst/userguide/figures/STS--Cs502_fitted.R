###################################################
### code chunk number 22: Cs502_fitted
###################################################
fitted2 <- tsSmooth(fit2, type = "xtt")
fitted2 <- subset(fitted2, .rownames %in% c("X1", "X2", "X3"))


