###################################################
### code chunk number 37: Cs_034_sleepmod4
###################################################
inits <- list(A = coef(sleep.mod3)$A, D = coef(sleep.mod3)$D)
# estimate a separate intercept for each but slope is the same
sleep.model <- list(
  A = "unequal", B = "diagonal and unequal", x0 = "zero", U = "zero",
  D = "unequal", d = exp.var, tinitx = 0, Q = "diagonal and unequal",
  R = "diagonal and unequal"
)
sleep.mod4 <- MARSS(dat, model = sleep.model, inits = inits, silent = TRUE)


