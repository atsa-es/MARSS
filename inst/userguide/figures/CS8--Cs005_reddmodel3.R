###################################################
### code chunk number 7: Cs005_reddmodel3
###################################################
model3 <- list()
model3$Q <- "diagonal and equal"
model3$R <- "diagonal and equal"
model3$U <- "equal"
model3$Z <- "identity"
model3$A <- "zero"
kem3 <- MARSS(logRedds, model = model3)


