###################################################
### code chunk number 5: Cs03_data3
###################################################
turtlename <- "BigMama"
theTurtle <- which(loggerheadNoisy$turtle == turtlename)
dat <- loggerheadNoisy[theTurtle, 5:6]
dat <- t(dat) # transpose


