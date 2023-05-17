###################################################
### code chunk number 35: Cs_032_sleepmod2-lm
###################################################
sleep.lm2 <- lm(Reaction ~ 0 + Subject + Days:Subject, data = sleepstudy)


