###################################################
### code chunk number 29: Cs011_fit-bird-model-2
###################################################
model.b2 <- list()
model.b2$Q <- "diagonal and equal"
model.b2$R <- "diagonal and equal"
model.b2$Z <- "identity"
model.b2$A <- "zero"
model.b2$U <- "equal"
kem.b2 <- MARSS(birddat, model = model.b2)


