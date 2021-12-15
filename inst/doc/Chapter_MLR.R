###################################################
### code chunk number 2: Cs_000_required_libraries
###################################################
library(MARSS)
library(xtable)
library(lattice)
library(nlme)
library(stringr)
library(lme4)


###################################################
### code chunk number 3: Cs_001_example1_plot
###################################################
data(longley)
plot(longley$Year, longley$Employed, type = "l", ylab = "Employed", xlab = "")


###################################################
### code chunk number 4: Cs_002_example1-data
###################################################
data(longley)
Employed <- matrix(longley$Employed, nrow = 1)


###################################################
### code chunk number 5: Cs_003_example1-params
###################################################
longley.model <- list()


###################################################
### code chunk number 6: Cs_004_example1-params1
###################################################
longley.model$U <- longley.model$Q <- "zero"
longley.model$C <- "zero"
longley.model$B <- longley.model$Z <- "identity"
longley.model$x0 <- "zero"
longley.model$tinitx <- 0


###################################################
### code chunk number 7: Cs_005_example1-paramsR
###################################################
longley.model$R <- matrix("r")


###################################################
### code chunk number 8: Cs_006_example1-paramsD
###################################################
longley.model$A <- matrix("intercept")
longley.model$D <- matrix(c("GNP", "Pop"), nrow = 1)


###################################################
### code chunk number 9: Cs_007_example1-eVar
###################################################
longley.model$d <- rbind(longley$GNP, longley$Population)


###################################################
### code chunk number 10: Cs_008_example1-marss
###################################################
mod1 <- MARSS(Employed, model = longley.model)


###################################################
### code chunk number 11: Cs_009_example1-marss
###################################################
coef(mod1, type = "vector")


###################################################
### code chunk number 12: Cs_010_example1-lm
###################################################
mod1.lm <- lm(Employed ~ GNP + Population, data = longley)
coef(mod1.lm)


###################################################
### code chunk number 13: Cs_011_example2-params
###################################################
longley.ar1 <- longley.model
longley.ar1$B <- matrix("b")
longley.ar1$Q <- matrix("q")


###################################################
### code chunk number 14: Cs_012_example2-marss
###################################################
mod2 <- MARSS(Employed, model = longley.ar1)


###################################################
### code chunk number 15: Cs_013_example2-marss-with-inits
###################################################
inits <- list(A = coef(mod1)$A, D = coef(mod1)$D)
mod2 <- MARSS(Employed,
  model = longley.ar1, inits = inits,
  control = list(maxit = 1000)
)
ests.marss <- c(
  b = coef(mod2)$B, alpha = coef(mod2)$A,
  GNP = coef(mod2)$D[1], Population = coef(mod2)$D[2],
  logLik = logLik(mod2)
)


###################################################
### code chunk number 16: Cs_014_example2-gls
###################################################
library(nlme)
mod2.gls <- gls(Employed ~ GNP + Population,
  correlation = corAR1(), data = longley, method = "ML"
)
mod2.gls.phi <- coef(mod2.gls$modelStruct[[1]], unconstrained = FALSE)
ests.gls <- c(
  b = mod2.gls.phi, alpha = coef(mod2.gls)[1],
  GNP = coef(mod2.gls)[2], Population = coef(mod2.gls)[3],
  logLik = logLik(mod2.gls)
)


###################################################
### code chunk number 17: Cs_014b_compare_marss_gls
###################################################
rbind(MARSS = ests.marss, GLS = ests.gls)


###################################################
### code chunk number 18: Cs_015_example2-plot
###################################################
pairs(longley)


###################################################
### code chunk number 19: Cs_016_full-model-list
###################################################
eVar.names <- colnames(longley)[-7]
eVar <- t(longley[, eVar.names])
longley.model <- list()
longley.model$U <- longley.model$Q <- "zero"
longley.model$C <- "zero"
longley.model$B <- longley.model$Z <- "identity"
longley.model$A <- matrix("intercept")
longley.model$R <- matrix("r")
longley.model$D <- matrix(eVar.names, nrow = 1)
longley.model$d <- eVar
longley.model$x0 <- "zero"
longley.model$tinitx <- 0


###################################################
### code chunk number 20: Cs_017_full-model-fit
###################################################
mod3.em <- MARSS(Employed, model = longley.model)
mod3.bfgs <- MARSS(Employed, model = longley.model, method = "BFGS")


###################################################
### code chunk number 21: Cs_018_full-em-ests
###################################################
par.names <- c("A.intercept", paste("D", eVar.names, sep = "."))
c(coef(mod3.em, type = "vector")[par.names], logLik = mod3.em$logLik)


###################################################
### code chunk number 22: Cs_019_full-bfgs-ests
###################################################
c(coef(mod3.bfgs, type = "vector")[par.names], logLik = mod3.bfgs$logLik)


###################################################
### code chunk number 23: Cs_020_full-lm-ests
###################################################
mod3.lm <- lm(Employed ~ 1 + GNP.deflator + GNP + Unemployed
  + Armed.Forces + Population + Year, data = longley)
c(coef(mod3.lm), logLik = logLik(mod3.lm))


###################################################
### code chunk number 24: Cs_021_full-model-correrr
###################################################
longley.correrr.model <- longley.model
longley.correrr.model$B <- matrix("b")
longley.correrr.model$Q <- matrix("q")


###################################################
### code chunk number 25: Cs_022_full-correrr-fit-hide
###################################################
inits <- list(A = coef(mod3.em)$A, D = coef(mod3.em)$D)
mod4.em <- MARSS(Employed, model = longley.correrr.model, inits = inits)
mod4.bfgs <- MARSS(Employed,
  model = longley.correrr.model,
  inits = inits, method = "BFGS"
)


###################################################
### code chunk number 26: Cs_023_full-correrr-em-ests
###################################################
c(coef(mod4.em, type = "vector")["B.b"], logLik = mod4.em$logLik)


###################################################
### code chunk number 27: Cs_024_full-correrr-bfgs-ests
###################################################
c(coef(mod4.bfgs, type = "vector")["B.b"], logLik = mod4.bfgs$logLik)


###################################################
### code chunk number 28: Cs_025_full-gls-ests
###################################################
mod4.gls <- gls(Employed ~ 1 + GNP.deflator + GNP + Unemployed
  + Armed.Forces + Population + Year,
correlation = corAR1(), data = longley, method = "ML"
)
mod4.gls.phi <- coef(mod4.gls$modelStruct[[1]], unconstrained = FALSE)
c(mod4.gls.phi, logLik = logLik(mod4.gls))


###################################################
### code chunk number 29: Cs_026_loadsleep
###################################################
data(sleepstudy, package = "lme4")


###################################################
### code chunk number 30: Cs_027_sleep-plot
###################################################
library(lattice)
xyplot(Reaction ~ Days | Subject, sleepstudy,
  type = c("g", "p", "r"),
  index = function(x, y) coef(lm(y ~ x))[1],
  xlab = "Days of sleep deprivation",
  ylab = "Average reaction time (ms)", aspect = "xy"
)


###################################################
### code chunk number 31: Cs_028_setupdata
###################################################
# number of subjects
nsub <- length(unique(sleepstudy$Subject))
ndays <- length(sleepstudy$Days) / nsub
dat <- matrix(sleepstudy$Reaction, nsub, ndays, byrow = TRUE)
rownames(dat) <- paste("sub", unique(sleepstudy$Subject), sep = ".")
exp.var <- matrix(sleepstudy$Days, 1, ndays, byrow = TRUE)


###################################################
### code chunk number 32: Cs_029_sleepmod1
###################################################
sleep.model <- list(
  A = "unequal", B = "zero", x0 = "zero", U = "zero",
  D = matrix("b1", nsub, 1), d = exp.var, tinitx = 0, Q = "zero"
)
sleep.mod1 <- MARSS(dat, model = sleep.model)


###################################################
### code chunk number 33: Cs_030_sleepmod1-lm
###################################################
sleep.lm1 <- lm(Reaction ~ -1 + Subject + Days, data = sleepstudy)


###################################################
### code chunk number 34: Cs_031_sleepmod2
###################################################
sleep.model <- list(
  A = "unequal", B = "zero", x0 = "zero", U = "zero",
  D = "unequal", d = exp.var, tinitx = 0, Q = "zero"
)
sleep.mod2 <- MARSS(dat, model = sleep.model, silent = TRUE)


###################################################
### code chunk number 35: Cs_032_sleepmod2-lm
###################################################
sleep.lm2 <- lm(Reaction ~ 0 + Subject + Days:Subject, data = sleepstudy)


###################################################
### code chunk number 36: Cs_033_sleepmod3
###################################################
sleep.model <- list(
  A = "unequal", B = "zero", x0 = "zero", U = "zero",
  D = "unequal", d = exp.var, tinitx = 0, Q = "zero",
  R = "diagonal and unequal"
)
sleep.mod3 <- MARSS(dat, model = sleep.model, silent = TRUE)


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


###################################################
### code chunk number 39: Cs_036_mod5-gls
###################################################
sleep.mod5.gls <- gls(Reaction ~ 0 + Subject + Days:Subject,
  data = sleepstudy,
  correlation = corAR1(form = ~ 1 | Subject), method = "ML"
)


###################################################
### code chunk number 40: Cs_037_makemodeltable
###################################################
if (!exists("tabledir")) tabledir <- ""
slope.names <- paste("D", rownames(dat), sep = ".")
phi.names <- names(coef(sleep.mod4, type = "vector"))[str_detect(names(coef(sleep.mod4, type = "vector")), "B.")]

model.data <- cbind(
  c(logLik(sleep.lm2), coef(sleep.lm2)[19:36], rep(NA, nsub)),
  c(sleep.mod2$logLik, coef(sleep.mod2, type = "vector")[slope.names], rep(NA, nsub)),
  c(sleep.mod3$logLik, coef(sleep.mod3, type = "vector")[slope.names], rep(NA, nsub)),
  c(sleep.mod4$logLik, coef(sleep.mod4, type = "vector")[c(slope.names, phi.names)]),
  c(sleep.mod5$logLik, coef(sleep.mod5, type = "vector")[c(slope.names, rep("B.diag", nsub))]),
  c(logLik(sleep.mod5.gls), coef(sleep.mod5.gls)[19:36], rep(coef(sleep.mod5.gls$modelStruct[[1]], unconstrained = FALSE), nsub))
)
rownames(model.data) <- c("logLik", paste("slope", unique(sleepstudy$Subject)), paste("phi", unique(sleepstudy$Subject)))
colnames(model.data) <- c("lm", "mod2 em", "mod3 em", "mod4 em", "mod5 em", "mod5 gls")
tmpaln <- "c" # figure out the number of cols automatically
for (i in 1:ncol(model.data)) tmpaln <- paste(tmpaln, "c", sep = "")
thetable <- xtable(model.data, caption = "Parameter estimates of different versions of the model where each subject has a separate intercept (response time on normal sleep) and different slope by day (increase in response time with each day of sleep deprivation).  The model types are discussed in the text.", label = "ref:tablesleepstudy", align = tmpaln, digits = 2)
print(thetable, type = "latex", file = paste(tabledir, "tablesleepstudy.tex", sep = ""), include.rownames = TRUE, include.colnames = TRUE, caption.placement = "top", table.placement = "htp", sanitize.text.function = function(x) {
  x
}, hline.after = c(-1, 0, nrow(model.data)))


