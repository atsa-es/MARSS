###################################################
### code chunk number 36: Cs_033_sleepmod3
###################################################
sleep.model <- list(
  A = "unequal", B = "zero", x0 = "zero", U = "zero",
  D = "unequal", d = exp.var, tinitx = 0, Q = "zero",
  R = "diagonal and unequal"
)
sleep.mod3 <- MARSS(dat, model = sleep.model, silent = TRUE)


