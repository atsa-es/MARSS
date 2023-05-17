###################################################
### code chunk number 34: Cs_031_sleepmod2
###################################################
sleep.model <- list(
  A = "unequal", B = "zero", x0 = "zero", U = "zero",
  D = "unequal", d = exp.var, tinitx = 0, Q = "zero"
)
sleep.mod2 <- MARSS(dat, model = sleep.model, silent = TRUE)


