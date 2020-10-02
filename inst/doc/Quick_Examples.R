###################################################
### code chunk number 2: Cs00_load_library
###################################################
library(MARSS)


###################################################
### code chunk number 3: Cs01_model-gen-spec
###################################################
Z <- matrix(list("z1", "z2", 0, 0, "z2", 3), 3, 2)
A <- matrix(0, 3, 1)
R <- matrix(list(0), 3, 3)
diag(R) <- c("r", "r", 1)
B <- matrix(list("b1", 0.1, "b2", 2), 2, 2)
U <- matrix(c("u", "u"), 2, 1)
Q <- matrix(c("q1", "q3", "q3", "q2"), 2, 2)
x0 <- matrix(c("pi1", "pi2"), 2, 1)
V0 <- diag(1, 2)
model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, 
                  Q = Q, x0 = x0, V0 = V0, tinitx = 0)


###################################################
### code chunk number 4: Cs02_model-general (eval = FALSE)
###################################################
## kemfit <- MARSS(dat, model = model.gen)


###################################################
### code chunk number 5: Cs03_enterdata
###################################################
dat <- t(harborSealWA)
dat <- dat[2:nrow(dat), ] # remove the year row


###################################################
### code chunk number 6: Cs04_model-default
###################################################
kemfit <- MARSS(dat)


###################################################
### code chunk number 7: Cs05_B-setup (eval = FALSE)
###################################################
## B <- Z <- diag(1, 5)
## U <- matrix(c("u1", "u2", "u3", "u4", "u5"), 5, 1)
## x0 <- A <- matrix(0, 5, 1)
## R <- Q <- matrix(list(0), 5, 5)
## diag(R) <- "r"
## diag(Q) <- c("q1", "q2", "q3", "q4", "q5")


###################################################
### code chunk number 8: Cs06_model-default-time
###################################################
kemfit.time <- system.time(MARSS(dat, silent = TRUE))


###################################################
### code chunk number 9: Cs07_model-bfgs
###################################################
kemfit.bfgs <- MARSS(dat, method = "BFGS")


###################################################
### code chunk number 10: Cs08_model-bfgs-time
###################################################
bfgsfit.time <- system.time(MARSS(dat, silent = TRUE, method = "BFGS"))


###################################################
### code chunk number 11: Cs09_model-bfgs2
###################################################
kemfit.bfgs2 <- MARSS(dat, method = "BFGS", inits = kemfit$par)


###################################################
### code chunk number 12: Cs10_model-default (eval = FALSE)
###################################################
## kemfit.strange <- MARSS(dat, model = list(tinitx = 1))


###################################################
### code chunk number 13: Cs11_model-corr1
###################################################
kemfit <- MARSS(dat, model = list(Q = "unconstrained"))


###################################################
### code chunk number 14: Cs12_model-corr2
###################################################
kemfit <- MARSS(dat, model = list(Q = "equalvarcov"))


###################################################
### code chunk number 15: Cs13_model-u-NS
###################################################
regions <- list("N", "N", "N", "S", "S")
U <- matrix(regions, 5, 1)
Q <- matrix(list(0), 5, 5)
diag(Q) <- regions
kemfit <- MARSS(dat, model = list(U = U, Q = Q))


###################################################
### code chunk number 16: Cs14_model-u-NS-fixR
###################################################
regions <- list("N", "N", "N", "S", "S")
U <- matrix(regions, 5, 1)
Q <- matrix(list(0), 5, 5)
diag(Q) <- regions
R <- diag(0.01, 5)
kemfit <- MARSS(dat, model = list(U = U, Q = Q, R = R))


###################################################
### code chunk number 17: Cs15_model-pan1
###################################################
Z <- factor(c(1, 1, 1, 1, 1))
kemfit <- MARSS(dat, model = list(Z = Z))


###################################################
### code chunk number 18: Cs16_model-pan2
###################################################
Z <- factor(c(1, 1, 1, 1, 1))
R <- "diagonal and unequal"
kemfit <- MARSS(dat, model = list(Z = Z, R = R))


###################################################
### code chunk number 19: Cs17_model-two1
###################################################
Z <- factor(c("N", "N", "N", "S", "S"))
Q <- "diagonal and equal"
U <- "equal"
kemfit <- MARSS(dat, model = list(Z = Z, Q = Q, U = U))


###################################################
### code chunk number 21: Cs17_model-two2
###################################################
kemfit <- MARSS(dat, model = list(Z = Z, Q = Q, U = U, A="zero"))


###################################################
### code chunk number 22: Cs18_model-time-varying
###################################################
U1 <- matrix("t1", 5, 1)
U2 <- matrix("t2", 5, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut[, , 1:floor(TT / 2)] <- U1
kemfit.tv <- MARSS(dat, model = list(U = Ut, Q = "diagonal and equal"))


###################################################
### code chunk number 24: Cs19_model-print
###################################################
print(kemfit)
print(kemfit$model)


###################################################
### code chunk number 23: Cs19b_model-time-varying2
###################################################
U1 <- matrix(c(rep("t1", 4), "hc"), 5, 1)
U2 <- matrix(c(rep("t2", 4), "hc"), 5, 1)
Ut <- array(U2, dim = c(dim(U1), dim(dat)[2]))
Ut[, , 1:floor(TT / 2)] <- U1
kemfit.tv <- MARSS(dat, model = list(U = Ut, Q = "diagonal and equal"))


###################################################
### code chunk number 25: Cs20_model-summary
###################################################
summary(kemfit$model)


###################################################
### code chunk number 26: Cs21_model-print-par
###################################################
print(kemfit, what = "par")


###################################################
### code chunk number 27: Cs22_model-print-Q
###################################################
print(kemfit, what = "Q")


###################################################
### code chunk number 28: Cs23_model-print-R
###################################################
x <- print(kemfit, what = "states", silent = TRUE)


###################################################
### code chunk number 29: Cs24_model-tidy-R
###################################################
broom::tidy(kemfit)
broom::glance(kemfit)


###################################################
### code chunk number 30: Cs25_CIs-hessian
###################################################
kem.with.hess.CIs <- MARSSparamCIs(kemfit)


###################################################
### code chunk number 31: Cs26_print-CIs
###################################################
print(kem.with.hess.CIs)


###################################################
### code chunk number 32: Cs27_CIs-pboot
###################################################
kem.w.boot.CIs <- MARSSparamCIs(kemfit, method = "parametric", nboot = 10)
# nboot should be more like 1000, but set low for example's sake
print(kem.w.boot.CIs)


###################################################
### code chunk number 33: Cs28_parvec
###################################################
parvec <- coef(kemfit, type = "vector")
parvec


###################################################
### code chunk number 34: Cs29_tsSmooth
###################################################
df <- tsSmooth(kemfit)
head(df)


###################################################
### code chunk number 35: Cs30_marsskf
###################################################
kf <- MARSSkf(kemfit)
names(kf)
# if you only need the logLik, 
MARSSkf(kemfit, only.logLik = TRUE)
# or
logLik(kemfit)


###################################################
### code chunk number 36: Cs31_like.kem.degen
###################################################
dat.short <- dat[1:4, 1:10]
kem.degen <- MARSS(dat.short, control = list(allow.degen = FALSE))


###################################################
### code chunk number 37: Cs32_like.kem200 (eval = FALSE)
###################################################
## kem.degen2 <- MARSS(dat.short, control = list(
##   maxit = 1000,
##   allow.degen = FALSE
## ), silent = 2)


###################################################
### code chunk number 38: Cs33_like.kem.allow.degen
###################################################
kem.short <- MARSS(dat.short)


###################################################
### code chunk number 39: Cs34_like.kem200
###################################################
kem.small <- MARSS(dat.short, model = list(
  Q = "diagonal and equal",
  U = "equal"
))


###################################################
### code chunk number 40: Cs35_boot
###################################################
boot.params <- MARSSboot(kemfit,
  nboot = 20, output = "parameters", sim = "parametric"
)$boot.params


###################################################
### code chunk number 41: Cs36_sim
###################################################
sim.data <- MARSSsimulate(kemfit, nsim = 2, tSteps = 100)$sim.data


###################################################
### code chunk number 42: Cs37_sim-fit
###################################################
kem.sim.1 <- MARSS(sim.data[, , 1])


###################################################
### code chunk number 43: Cs38_sim-like
###################################################
kem.sim.2 <- kem.sim.1
kem.sim.2$model$data <- sim.data[, , 2]
MARSSkf(kem.sim.2)$logLik


###################################################
### code chunk number 44: Cs39_AICb
###################################################
kemfit.with.AICb <- MARSSaic(kemfit,
  output = "AICbp",
  Options = list(nboot = 10, silent = TRUE)
)

print(kemfit.with.AICb)


