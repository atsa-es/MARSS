###################################################
### code chunk number 2: Cs_01_readindata
###################################################
data(SalmonSurvCUI)
years <- SalmonSurvCUI[, 1]
TT <- length(years)
# response data: logit(survival)
dat <- matrix(SalmonSurvCUI[, 2], nrow = 1)


###################################################
### code chunk number 3: Cs_02_zscore
###################################################
CUI <- SalmonSurvCUI[, "CUI.apr"]
CUI.z <- zscore(CUI)
# number of state = # of regression params (slope(s) + intercept)
m <- 1 + 1


###################################################
### code chunk number 4: Cs_030_plotdata
###################################################
par(mfrow = c(m, 1), mar = c(4, 4, 0.1, 0), oma = c(0, 0, 2, 0.5))
plot(years, dat, xlab = "", ylab = "Logit(s)", bty = "n", xaxt = "n", pch = 16, col = "darkgreen", type = "b")
plot(years, CUI.z, xlab = "", ylab = "CUI", bty = "n", xaxt = "n", pch = 16, col = "blue", type = "b")
axis(1, at = seq(1965, 2005, 5))
mtext("Year of ocean entry", 1, line = 3)


###################################################
### code chunk number 5: Cs_031_univDLMproc
###################################################
# for process eqn
B <- diag(m) # 2x2; Identity
U <- matrix(0, nrow = m, ncol = 1) # 2x1; both elements = 0
Q <- matrix(list(0), m, m) # 2x2; all 0 for now
diag(Q) <- c("q1", "q2") # 2x2; diag = (q1,q2)


###################################################
### code chunk number 6: Cs_04_univDLMobs
###################################################
# for observation eqn
Z <- array(NA, c(1, m, TT)) # NxMxT; empty for now
Z[1, 1, ] <- rep(1, TT) # Nx1; 1's for intercept
Z[1, 2, ] <- CUI.z # Nx1; regr variable
A <- matrix(0) # 1x1; scalar = 0
R <- matrix("r") # 1x1; scalar = r


###################################################
### code chunk number 7: Cs_05_univDLM-list
###################################################
# only need starting values for regr parameters
inits.list <- list(x0 = matrix(c(0, 0), nrow = m))
# list of model matrices & vectors
mod.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)


###################################################
### code chunk number 8: Cs_06_univDLM-fit
###################################################
dlm1 <- MARSS(dat, inits = inits.list, model = mod.list)


###################################################
### code chunk number 9: Cs_07_plotdlm1
###################################################
ylabs <- c(expression(alpha[t]), expression(beta[t]))
colr <- c("darkgreen", "blue")
par(mfrow = c(m, 1), mar = c(4, 4, 0.1, 0), oma = c(0, 0, 2, 0.5))
for (i in 1:m) {
  mn <- dlm1$states[i, ]
  se <- dlm1$states.se[i, ]
  plot(years, mn,
    xlab = "", ylab = ylabs[i], bty = "n", xaxt = "n", type = "n",
    ylim = c(min(mn - 2 * se), max(mn + 2 * se))
  )
  lines(years, rep(0, TT), lty = "dashed")
  lines(years, mn, col = colr[i], lwd = 3)
  lines(years, mn + 2 * se, col = colr[i])
  lines(years, mn - 2 * se, col = colr[i])
}
axis(1, at = seq(1965, 2005, 5))
mtext("Year of ocean entry", 1, line = 3)


###################################################
### code chunk number 10: Cs_08_univDLM-fore-mean
###################################################
# get list of Kalman filter output
kf.out <- MARSSkfss(dlm1)
# forecasts of regr parameters; 2xT matrix
eta <- kf.out$xtt1
# ts of E(forecasts)
fore.mean <- vector()
for (t in 1:TT) {
  fore.mean[t] <- Z[, , t] %*% eta[, t, drop = F]
}


###################################################
### code chunk number 11: Cs_09_univDLM-fore-Var
###################################################
# variance of regr parameters; 1x2xT array
Phi <- kf.out$Vtt1
# obs variance; 1x1 matrix
R.est <- coef(dlm1, type = "matrix")$R
# ts of Var(forecasts)
fore.var <- vector()
for (t in 1:TT) {
  tZ <- matrix(Z[, , t], m, 1) # transpose of Z
  fore.var[t] <- Z[, , t] %*% Phi[, , t] %*% tZ + R.est
}


###################################################
### code chunk number 12: Cs_10_plot-dlm-Forecast-Logit
###################################################
par(mar = c(4, 4, 0.1, 0), oma = c(0, 0, 2, 0.5))
ylims <- c(min(fore.mean - 2 * sqrt(fore.var)), max(fore.mean + 2 * sqrt(fore.var)))
plot(years, t(dat),
  type = "p", pch = 16, ylim = ylims,
  col = "blue", xlab = "", ylab = "Logit(s)", xaxt = "n"
)
lines(years, fore.mean, type = "l", xaxt = "n", ylab = "", lwd = 3)
lines(years, fore.mean + 2 * sqrt(fore.var))
lines(years, fore.mean - 2 * sqrt(fore.var))
axis(1, at = seq(1965, 2005, 5))
mtext("Year of ocean entry", 1, line = 3)


###################################################
### code chunk number 13: Cs_11_plot-dlm-Forecast-Raw
###################################################
invLogit <- function(x) {
  1 / (1 + exp(-x))
}
ff <- invLogit(fore.mean)
fup <- invLogit(fore.mean + 2 * sqrt(fore.var))
flo <- invLogit(fore.mean - 2 * sqrt(fore.var))
par(mar = c(4, 4, 0.1, 0), oma = c(0, 0, 2, 0.5))
ylims <- c(min(flo), max(fup))
plot(years, invLogit(t(dat)),
  type = "p", pch = 16, ylim = ylims,
  col = "blue", xlab = "", ylab = "Survival", xaxt = "n"
)
lines(years, ff, type = "l", xaxt = "n", ylab = "", lwd = 3)
lines(years, fup)
lines(years, flo)
axis(1, at = seq(1965, 2005, 5))
mtext("Year of ocean entry", 1, line = 3)


###################################################
### code chunk number 14: Cs_12_dlm-forecast-errors
###################################################
# forecast errors
innov <- kf.out$Innov


###################################################
### code chunk number 16: Cs_13_plot-dlmQQ
###################################################
# use layout to get nicer plots
layout(matrix(c(0, 1, 1, 1, 0), 1, 5, byrow = TRUE))
# set up L plotting space
par(mar = c(4, 4, 1, 0), oma = c(0, 0, 0, 0.5))
# Q-Q plot of innovations
qqnorm(t(innov), main = "", pch = 16, col = "blue")
qqline(t(innov))
# set up R plotting space
# par(mar=c(4,0,1,1)) #, oma=c(0,0,0,0.5))
# boxplot of innovations
# boxplot(t(innov), axes=FALSE)


###################################################
### code chunk number 17: Cs_14_dlm-Innov-t-test
###################################################
# p-value for t-test of H0: E(innov) = 0
t.test(t(innov), mu = 0)$p.value


###################################################
### code chunk number 19: Cs_15_plot-dlm-ACF
###################################################
# use layout to get nicer plots
layout(matrix(c(0, 1, 1, 1, 0), 1, 5, byrow = TRUE))
# set up plotting space
par(mar = c(4, 4, 1, 0), oma = c(0, 0, 0, 0.5))
# ACF of innovations
acf(t(innov), lwd = 2, lag.max = 10)


