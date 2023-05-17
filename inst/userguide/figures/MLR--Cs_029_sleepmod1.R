###################################################
### code chunk number 32: Cs_029_sleepmod1
###################################################
sleep.model <- list(
  A = "unequal", B = "zero", x0 = "zero", U = "zero",
  D = matrix("b1", nsub, 1), d = exp.var, tinitx = 0, Q = "zero"
)
sleep.mod1 <- MARSS(dat, model = sleep.model)


