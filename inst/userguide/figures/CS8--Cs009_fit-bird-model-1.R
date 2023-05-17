###################################################
### code chunk number 28: Cs009_fit-bird-model-1
###################################################
model.b1=list()
model.b1$R="diagonal and equal"
model.b1$Z=matrix(1,3,1)
kem.b1 = MARSS(birddat, model=model.b1, control=list(minit=100) )


