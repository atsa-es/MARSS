###################################################
### code chunk number 10: Cs2_Code1
###################################################
# Code to fit the single population model with i.i.d. errors
# Read in data
dat <- t(harborSealWA) # MARSS needs time ACROSS columns
years <- dat[1, ]
n <- nrow(dat) - 1
dat <- dat[2:nrow(dat), ]
legendnames <- (unlist(dimnames(dat)[1]))

# estimate parameters
Z.model <- factor(c(1, 1, 1, 1, 1))
R.model <- "diagonal and equal"
kem1 <- MARSS(dat, model = list(Z = Z.model, R = R.model))

# make figure
graphics::matplot(years, t(dat),
  xlab = "", ylab = "index of log abundance",
  pch = c("1", "2", "3", "4", "5"), ylim = c(5, 9), bty = "L"
)
lines(years, kem1$states - 1.96 * kem1$states.se,
  type = "l",
  lwd = 1, lty = 2, col = "red"
)
lines(years, kem1$states + 1.96 * kem1$states.se,
  type = "l",
  lwd = 1, lty = 2, col = "red"
)
lines(years, kem1$states, type = "l", lwd = 2)
title("Observations and total population estimate", cex.main = .9)

coef(kem1, type = "vector") # parameter estimates as a vector

# show estimated elements for each parameter matrix as a list
coef(kem1)

kem1$logLik # show the log-likelihood
kem1$AIC # show the AIC


###################################################
### code chunk number 17: Cs2_Code2
###################################################
# Fit the single population model with independent and unequal errors
Z.model <- factor(c(1, 1, 1, 1, 1))
R.model <- "diagonal and unequal"
kem2 <- MARSS(dat, model = list(Z = Z.model, R = R.model))

coef(kem2) # the estimated parameter elements
kem2$logLik # log likelihood
c(kem1$AIC, kem2$AIC) # AICs

resids <- MARSSresiduals(kem2, type = "tt1")$model.residuals
par(mfrow = c(2, 3))
for (i in 1:n) {
  plot(resids[i,], ylab = "residuals")
  title(legendnames[i])
}
par(mfrow = c(1, 1))


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
# plot residuals
resids <- MARSSresiduals(kem3, type = "tt1")$model.residuals
par(mfrow = c(2, 3))
for (i in 1:n) {
  plot(resids[i,], ylab = "residuals")
  title(legendnames[i])
}
par(mfrow = c(1, 1))


###################################################
### code chunk number 22: Cs2_Code4
###################################################
Z.model <- factor(c(1, 2, 3, 4, 5))
U.model <- "equal"
Q.model <- "diagonal and equal"
R.model <- "diagonal and unequal"
kem <- MARSS(dat, model = list(
  Z = Z.model,
  U = U.model, Q = Q.model, R = R.model
))


###################################################
### code chunk number 28: Cs2_Code5_7
###################################################
# Two subpopulations with different population parameters
Z.model <- factor(c(1, 1, 2, 2, 2))
U.model <- "unequal"
Q.model <- "diagonal and unequal"
R.model <- "diagonal and unequal"
kem <- MARSS(dat, model = list(Z = Z.model, U = U.model, Q = Q.model, R = R.model))

# Hood Canal covaries with the other regions
Z.model <- factor(c(1, 1, 1, 1, 2))
U.model <- "unequal"
Q.model <- "equalvarcov"
R.model <- "diagonal and unequal"
kem <- MARSS(dat, model = list(Z = Z.model, U = U.model, Q = Q.model, R = R.model))

# Three subpopulations with shared parameter values
Z.model <- factor(c(1, 1, 2, 2, 3))
U.model <- "unequal"
Q.model <- "diagonal and unequal"
R.model <- "diagonal and unequal"
kem <- MARSS(dat, model = list(Z = Z.model, U = U.model, Q = Q.model, R = R.model))


