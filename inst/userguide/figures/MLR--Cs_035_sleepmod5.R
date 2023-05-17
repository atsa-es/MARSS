###################################################
### code chunk number 38: Cs_035_sleepmod5
###################################################
inits <- list(A = coef(sleep.mod3)$A, D = coef(sleep.mod3)$D)
# estimate a separate intercept for each but slope is the same
sleep.model <- list(
  A = "unequal", B = "diagonal and equal", x0 = "zero", U = "zero",
  D = "unequal", d = exp.var, tinitx = 0, Q = "diagonal and equal",
  R = "diagonal and equal"
)
sleep.mod5 <- MARSS(dat, model = sleep.model, inits = inits, silent = TRUE)


