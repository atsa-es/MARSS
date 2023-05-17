###################################################
### code chunk number 20: Cs2_Code3
###################################################
# fit the north and south population model
Z.model <- factor(c(1, 1, 2, 2, 2))
U.model <- "equal"
Q.model <- "diagonal and equal"
R.model <- "diagonal and unequal"
kem3 <- MARSS(dat, model = list(
  Z = Z.model,
  R = R.model, U = U.model, Q = Q.model
))
# plot smoothation residuals
plot(kem3, plot.type="model.resids.ytT")


