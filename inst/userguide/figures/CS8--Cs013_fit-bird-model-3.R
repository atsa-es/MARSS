###################################################
### code chunk number 30: Cs013_fit-bird-model-3
###################################################
model.b3 <- model.b2 # is is based on model.b2
# all we change is the structure of Q
model.b3$Q <- "diagonal and unequal"
model.b3$U <- "unequal"
kem.b3 <- MARSS(birddat, model = model.b3)


