###################################################
### code chunk number 5: Cs003_reddmodel1
###################################################
model1 <- list()
model1$R <- "diagonal and equal"
model1$Z <- matrix(1, 2, 1)
model1$A <- "scaling"
kem1 <- MARSS(logRedds, model = model1)


