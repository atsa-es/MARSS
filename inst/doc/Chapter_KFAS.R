###################################################
### code chunk number 2: Cs00_required-libraries
###################################################
library(MARSS)
library(KFAS)
library(ggplot2) # plotting
library(tidyr) # data frame manipulation


###################################################
### code chunk number 3: Cs101_fitting-models
###################################################
model_Nile <- SSModel(Nile ~ SSMtrend(
  degree = 1,
  Q = list(matrix(NA))
),
H = matrix(NA)
)
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_kfas_default <- fitSSM(model_Nile, kinits, method = "BFGS")


###################################################
### code chunk number 4: Cs102_fitting-models
###################################################
model_Nile_stoch <- model_Nile
model_Nile_stoch$a1[1, 1] <- 0
model_Nile_stoch$P1[1, 1] <- 1000
model_Nile_stoch$P1inf[1, 1] <- 0
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_kfas_stoch <- fitSSM(model_Nile_stoch, kinits, method = "BFGS")
kfs_kfas_stoch <- KFS(fit_kfas_stoch$model)


###################################################
### code chunk number 5: Cs103_fitting-models
###################################################
mod.nile <- list(
  Z = matrix(1), A = matrix(0), R = matrix("r"),
  B = matrix(1), U = matrix(0), Q = matrix("q"),
  tinitx = 1
)


###################################################
### code chunk number 6: Cs104_fitting-models
###################################################
dat <- t(as.matrix(Nile))
rownames(dat) <- "Nile"
fit_em_default <- MARSS(dat, model = mod.nile, silent = TRUE)
inits <- list(Q = matrix(var(Nile)), R = matrix(var(Nile)))
fit_bfgs_default <- MARSS(dat,
  model = mod.nile, inits = inits,
  method = "BFGS", silent = TRUE
)


###################################################
### code chunk number 7: Cs105_fitting-models
###################################################
mod.nile.stoch <- mod.nile
mod.nile.stoch$x0 <- fit_kfas_stoch$model$a1
mod.nile.stoch$V0 <- fit_kfas_stoch$model$P1
fit_em_stoch <- MARSS(dat, model = mod.nile.stoch, silent = TRUE)
fit_bfgs_stoch <- MARSS(dat,
  model = mod.nile.stoch, inits = inits,
  method = "BFGS", silent = TRUE
)


###################################################
### code chunk number 8: Cs106_fitting-models
###################################################
marss_kfas_model <- MARSSkfas(fit_em_stoch,
  return.kfas.model = TRUE,
  return.lag.one = FALSE
)$kfas.model
marss_kfas_model$Q[1, 1, 1] <- NA
marss_kfas_model$H[1, 1, 1] <- NA
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_marss_kfas <- fitSSM(marss_kfas_model, kinits, method = "BFGS")


###################################################
### code chunk number 9: Cs107_fitting-models
###################################################
vals <- rbind(
  c(fit_kfas_default$model$Q, fit_kfas_default$model$H, -1 * fit_kfas_default$optim.out$value),
  c(coef(fit_em_default)$Q, coef(fit_em_default)$R, logLik(fit_em_default)),
  c(coef(fit_bfgs_default)$Q, coef(fit_bfgs_default)$R, logLik(fit_bfgs_default)),
  c(fit_kfas_stoch$model$Q, fit_kfas_stoch$model$H, -1 * fit_kfas_stoch$optim.out$value),
  c(coef(fit_em_stoch)$Q, coef(fit_em_stoch)$R, logLik(fit_em_stoch)),
  c(coef(fit_bfgs_stoch)$Q, coef(fit_bfgs_stoch)$R, logLik(fit_bfgs_stoch)),
  c(fit_marss_kfas$model$Q, fit_marss_kfas$model$H, -1 * fit_marss_kfas$optim.out$value)
)
rownames(vals) <- c(
  "KFAS default", "MARSS em default", "MARSS bfgs default",
  "KFAS stoch", "MARSS em stoch", "MARSS bfgs stoch", "KFAS w marss kfas model"
)
colnames(vals) <- c("Q", "R", "logLik")
vals


###################################################
### code chunk number 10: Cs201_state-filtering
###################################################
fit_kfas <- fit_kfas_stoch
fit_marss <- fit_em_stoch
fit_marss$par$Q[1, 1] <- fit_kfas$model$Q
fit_marss$par$R[1, 1] <- fit_kfas$model$H


###################################################
### code chunk number 11: Cs202_state-filtering
###################################################
kf_kfas <- KFS(fit_kfas$model,
  filtering = "state",
  smoothing = "state", simplify = FALSE
)


###################################################
### code chunk number 12: Cs203_state-filtering
###################################################
kf_marss <- MARSSkfss(fit_marss)


###################################################
### code chunk number 13: Cs204_state-filtering
###################################################
names(kf_kfas)
names(kf_marss)


###################################################
### code chunk number 14: Cs205_state-filtering
###################################################
n <- 5


###################################################
### code chunk number 15: Cs206_state-filtering
###################################################
cbind(
  a = kf_kfas$a[1:n], xtt1 = kf_marss$xtt1[1:n],
  att = kf_kfas$att[1:n], xtt = kf_marss$xtt[1:n]
)


###################################################
### code chunk number 16: Cs207_state-filtering
###################################################
cbind(kf_kfas$alphahat[1:n], kf_marss$xtT[1:n])


###################################################
### code chunk number 17: Cs208_state-filtering
###################################################
cbind(
  v = kf_kfas$v[1:n], Innov = kf_marss$Innov[1:n],
  F = kf_kfas$F[1:n], Sigma = kf_marss$Sigma[1:n]
)


###################################################
### code chunk number 18: Cs209_state-filtering
###################################################
cbind(
  P = kf_kfas$P[1:n], Vtt1 = kf_marss$Vtt1[1:n],
  Ptt = kf_kfas$Ptt[1:n], Vtt = kf_marss$Vtt[1:n]
)


###################################################
### code chunk number 19: Cs301_obs-filtering
###################################################
f_kfas <- KFS(fit_kfas$model,
  filtering = "signal",
  smoothing = "signal", simplify = FALSE
)


###################################################
### code chunk number 20: Cs302_obs-filtering
###################################################
f_marss <- MARSSkf(fit_marss)


###################################################
### code chunk number 21: Cs304_obs-filtering
###################################################
n <- 10


###################################################
### code chunk number 22: Cs305_obs-filtering
###################################################
ytt1_fit <- fitted(fit_marss, type = "ytt1")$.fitted
ytt1_hatyt <- MARSShatyt(fit_marss, only.kem = FALSE)$ytt1
cbind(m = f_kfas$m[1:n], fitted = ytt1_fit[1:n], MARSShatyt = ytt1_hatyt[1:n])


###################################################
### code chunk number 23: Cs306_obs-filtering
###################################################
var.Eytt1_fit <-
  fitted(fit_marss, type = "ytt1", interval = "confidence")$.se^2
var.Eytt1_hatyt <-
  MARSShatyt(fit_marss, only.kem = FALSE)$var.Eytt1
cbind(
  P_mu = f_kfas$P_mu[1:n], fitted = var.Eytt1_fit[1:n],
  MARSShatyt = var.Eytt1_hatyt[1:n]
)


###################################################
### code chunk number 24: Cs307_obs-filtering
###################################################
ytT_fit <- fitted(fit_marss, type = "ytT")$.fitted
ytT_hatyt <- MARSShatyt(fit_marss)$ytT
cbind(
  a = f_kfas$muhat[1:n], fitted = ytT_fit[1:n],
  MARSShatyt = ytT_hatyt[1:n], Nile = Nile[1:n]
)


###################################################
### code chunk number 25: Cs308_obs-filtering
###################################################
var.EytT_fit <-
  fitted(fit_marss, type = "ytT", interval = "confidence")$.se^2
var.EytT_hatyt <-
  MARSShatyt(fit_marss, only.kem = FALSE)$var.EytT
cbind(
  V_mu = f_kfas$V_mu[1:n], fitted = var.EytT_fit[1:n],
  MARSShatyt = var.EytT_hatyt[1:n]
)


###################################################
### code chunk number 26: Cs401_conf-int
###################################################
conf_kfas <- predict(fit_kfas$model,
  interval = "confidence",
  se.fit = TRUE
)
head(conf_kfas)


###################################################
### code chunk number 27: Cs402_conf-int
###################################################
conf_marss1 <- fitted(fit_marss, type = "ytT", interval = "confidence")
head(conf_marss1)


###################################################
### code chunk number 28: Cs403_conf-int
###################################################
conf_marss2 <- predict(fit_marss,
  type = "ytT",
  interval = "confidence", level = 0.95
)
head(conf_marss2$pred)


###################################################
### code chunk number 29: Cs404_conf-int
###################################################
pred_kfas <- predict(fit_kfas$model,
  interval = "prediction", se.fit = TRUE
)
head(pred_kfas)


###################################################
### code chunk number 30: Cs405_conf-int
###################################################
pred_marss1 <- fitted(fit_marss, type = "ytT", interval = "prediction")
head(pred_marss1)


###################################################
### code chunk number 31: Cs406_conf-int
###################################################
pred_marss2 <- predict(fit_marss,
  type = "ytT",
  interval = "prediction", level = 0.95
)


###################################################
### code chunk number 32: Cs407_conf-int
###################################################
conf_kfas_t1 <- predict(fit_kfas$model,
  interval = "confidence",
  se.fit = TRUE, filtered = TRUE
)
head(conf_kfas_t1)


###################################################
### code chunk number 33: Cs408_conf-int
###################################################
conf_marss1_t1 <- fitted(fit_marss, type = "ytt1", interval = "confidence")
head(conf_marss1_t1)


###################################################
### code chunk number 34: Cs409_conf-int
###################################################
conf_marss2_t1 <- predict(fit_marss,
  type = "ytt1",
  interval = "confidence", level = 0.95
)
head(conf_marss2_t1$pred)


###################################################
### code chunk number 35: Cs410_conf-int
###################################################
pred_kfas_t1 <- predict(fit_kfas$model,
  interval = "prediction",
  se.fit = TRUE, filtered = TRUE
)
head(pred_kfas_t1)


###################################################
### code chunk number 36: Cs411_conf-int
###################################################
pred_marss1_t1 <- fitted(fit_marss, type = "ytt1", interval = "prediction")
head(pred_marss1_t1)


###################################################
### code chunk number 37: Cs412_conf-int
###################################################
pred_marss2_t1 <- predict(fit_marss,
  type = "ytt1",
  interval = "prediction", level = 0.95
)


###################################################
### code chunk number 49: Cs501_plotting
###################################################
ts.plot(cbind(Nile, pred_kfas[, c("fit", "lwr", "upr")], conf_kfas[, c("lwr", "upr")]),
  col = c(1:2, 3, 3, 4, 4),
  ylab = "Predicted Annual flow", main = "River Nile"
)


###################################################
### code chunk number 38: Cs501_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- residuals(kfs, type = "recursive")
resid_marss <- residuals(fit_marss,
  type = "tt1",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 50: Cs502_plotting
###################################################
plot(fit_marss, plot.type = "model.ytT", pi.int = TRUE)


###################################################
### code chunk number 39: Cs502_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- rstandard(kfs,
  type = "recursive",
  standardization_type = "marginal"
)
resid_marss <- residuals(fit_marss,
  type = "tt1",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 51: Cs503_plotting
###################################################
require(ggplot2)
df <- cbind(conf_marss1, pred_marss1[, c(".lwr", ".upr")])
ggplot(df, aes(x = t, y = .fitted)) +
  geom_ribbon(aes(ymin = .lwr, ymax = .upr), fill = "grey") +
  geom_ribbon(aes(ymin = .conf.low, ymax = .conf.up), fill = "blue", alpha = 0.25) +
  geom_line(linetype = 2) +
  ylab("Predicted Annual Flow") +
  xlab("") +
  ggtitle("River Nile")


###################################################
### code chunk number 40: Cs503_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- rstandard(kfs,
  type = "recursive",
  standardization_type = "cholesky"
)
resid_marss <- residuals(fit_marss,
  type = "tt1",
  standardization = "Cholesky"
)
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 41: Cs504_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- residuals(kfs, type = "pearson")
resid_marss <- residuals(fit_marss, type = "tT")
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 42: Cs505_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- rstandard(kfs,
  type = "pearson",
  standardization_type = "marginal"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 43: Cs506_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- rstandard(kfs,
  type = "pearson",
  standardization_type = "cholesky"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "Cholesky"
)
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 44: Cs507_residuals
###################################################
kfs <- KFS(fit_kfas$model)
resid_kfas <- residuals(kfs, type = "response")
resid_marss <- residuals(fit_marss, type = "tT")
resid_marss <- subset(resid_marss, name == "model")
df <- cbind(
  MARSS = resid_marss$.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 45: Cs508_residuals
###################################################
kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
resid_kfas <- residuals(kfs, type = "state")
resid_marss <- residuals(fit_marss, type = "tT")
resid_marss <- subset(resid_marss, name == "state")
df <- cbind(
  MARSS = resid_marss$.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 46: Cs509_residuals
###################################################
kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
resid_kfas <- rstandard(kfs,
  type = "state",
  standardization_type = "marginal"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "state")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 47: Cs509_residuals
###################################################
kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
resid_kfas <- rstandard(kfs,
  type = "state",
  standardization_type = "cholesky"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "Block.Cholesky"
)
resid_marss <- subset(resid_marss, name == "state")
df <- cbind(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas)
)
head(df)


###################################################
### code chunk number 48: Cs510_residuals
###################################################
kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
test <- cbind(
  b = fit_kfas$model$Q[1, 1, 1] - kfs$V_eta[1, 1, ],
  a = MARSSresiduals(fit_marss, type = "tT")$var.residuals[2, 2, ]
)
test <- as.data.frame(test)
test$diff <- test$b - test$a
head(test)
tail(test)


###################################################
### code chunk number 52: Cs601_missing-values
###################################################
NileNA <- Nile
NileNA[c(21:40, 61:80)] <- NA
model_NileNA_stoch <-
  SSModel(NileNA ~ SSMtrend(
    degree = 1,
    Q = list(matrix(NA))
  ),
  H = matrix(NA)
  )
model_NileNA_stoch$a1[1, 1] <- 0
model_NileNA_stoch$P1[1, 1] <- model_Nile_stoch$P1[1, 1]
model_NileNA_stoch$P1inf[1, 1] <- 0
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_kfas_NA <- fitSSM(model_NileNA_stoch, kinits, method = "BFGS")
fit_marss_NA <- MARSS(as.vector(NileNA),
  model = mod.nile.stoch,
  inits = inits, method = "BFGS", silent = TRUE
)


###################################################
### code chunk number 53: Cs602_missing-values
###################################################
rbind(
  MARSS = c(
    Q = coef(fit_marss_NA, type = "matrix")$Q,
    R = coef(fit_marss_NA, type = "matrix")$R,
    logLik = logLik(fit_marss_NA)
  ),
  KFAS = c(
    Q = fit_kfas_NA$model$Q,
    R = fit_kfas_NA$model$H,
    logLik = -1 * fit_kfas_NA$optim.out$value
  )
)


###################################################
### code chunk number 55: Cs603_missing-values
###################################################
conf_kfas_NA <-
  predict(fit_kfas_NA$model, interval = "confidence", filtered = FALSE)
conf_marss_NA <-
  predict(fit_marss_NA, interval = "confidence", type = "ytT", level = 0.95)$pred


###################################################
### code chunk number 56: Cs604_marss-mult-fig-2
###################################################
require(ggplot2)
df1 <- as.data.frame(conf_kfas_NA)
df1$name <- "KFAS"
df2 <- conf_marss_NA[, c("estimate", "Lo 95", "Hi 95")]
df2$name <- "MARSS"
colnames(df2) <- colnames(df1)
df <- rbind(df1, df2)
df$t <- as.vector(time(NileNA))
df$y <- conf_marss_NA$y
ggplot(df, aes(x = t, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey") +
  geom_line() +
  ylab("Predicted Annual Flow") +
  xlab("") +
  ggtitle("River Nile with 95% CIs on estimate") +
  facet_wrap(~name)


###################################################
### code chunk number 57: Cs605_missing-values
###################################################
fitted_kfas_NA <- data.frame(
  smooth = as.vector(fitted(fit_kfas_NA$model)),
  one.step.ahead = as.vector(fitted(fit_kfas_NA$model, filtered = TRUE)),
  name = "KFAS"
)
fitted_marss_NA <- data.frame(
  smooth = fitted(fit_marss_NA, type = "ytT")$.fitted,
  one.step.ahead = fitted(fit_marss_NA, type = "ytt1")$.fitted,
  name = "MARSS"
)


###################################################
### code chunk number 58: Cs606_missing-values
###################################################
require(ggplot2)
require(tidyr)
df <- rbind(fitted_kfas_NA, fitted_marss_NA)
df$t <- as.vector(time(NileNA))
df$y <- conf_marss_NA$y
df <- tidyr::pivot_longer(df, c(smooth, one.step.ahead), names_to = "type", values_to = "value")
ggplot(df, aes(x = t, y = value, col = type)) +
  geom_point(aes(x = t, y = y), col = "blue", size = 0.5, na.rm = TRUE) +
  geom_line() +
  ylab("Predicted Annual Flow") +
  xlab("") +
  ggtitle("River Nile - smoothed versus filtered") +
  facet_wrap(~name, ncol = 1)


###################################################
### code chunk number 59: Cs701_globaltemp
###################################################
data("GlobalTemp")
ts.plot(GlobalTemp, col = 1:2, main = "Two ts for Global Temperature")


###################################################
### code chunk number 60: Cs702_globaltemp
###################################################
data("GlobalTemp")
model_temp <- SSModel(GlobalTemp ~ SSMtrend(1, Q = NA, type = "common"),
  H = matrix(NA, 2, 2)
)
kinits <- chol(cov(GlobalTemp))[c(1, 4, 3)]
kinits <- c(0.5 * log(0.1), log(kinits[1:2]), kinits[3])
kfas_temp_default <- fitSSM(model_temp, kinits, method = "BFGS")
model_temp_stoch <- model_temp
model_temp_stoch$a1[1, 1] <- 0
model_temp_stoch$P1[1, 1] <- 1000 * max(diag(var(GlobalTemp)))
model_temp_stoch$P1inf[1, 1] <- 0
kfas_temp_stoch <- fitSSM(model_temp_stoch, kinits, method = "BFGS")


###################################################
### code chunk number 61: Cs703_globaltemp
###################################################
mod.list <- list(
  Z = matrix(1, 2, 1),
  R = matrix(c("r1", "c", "c", "r2"), 2, 2),
  U = matrix(0),
  A = matrix(0, 2, 1),
  tinitx = 1
)
marss_temp_default <- MARSS(t(GlobalTemp), model = mod.list)
mod.list$x0 <- kfas_temp_stoch$model$a1
mod.list$V0 <- kfas_temp_stoch$model$P1
marss_temp_stoch_em <- MARSS(t(GlobalTemp), model = mod.list)
# use inits from a short run of EM algorithm
inits <- MARSS(t(GlobalTemp),
  model = mod.list, control = list(maxit = 20),
  silent = TRUE
)
marss_temp_stoch_bfgs <- MARSS(t(GlobalTemp),
  model = mod.list,
  inits = inits, method = "BFGS"
)


###################################################
### code chunk number 62: Cs704_globaltemp
###################################################
vals <- rbind(
  c(kfas_temp_default$model$Q, kfas_temp_default$model$H[c(1, 2, 4)], -1 * kfas_temp_default$optim.out$value),
  c(coef(marss_temp_default)$Q, coef(marss_temp_default)$R, logLik(marss_temp_default)),
  c(kfas_temp_stoch$model$Q, kfas_temp_stoch$model$H[c(1, 2, 4)], -1 * kfas_temp_stoch$optim.out$value),
  c(coef(marss_temp_stoch_em)$Q, coef(marss_temp_stoch_em)$R, logLik(marss_temp_stoch_em)),
  c(coef(marss_temp_stoch_bfgs)$Q, coef(marss_temp_stoch_bfgs)$R, logLik(marss_temp_stoch_bfgs))
)
rownames(vals) <- c(
  "KFAS default", "MARSS em default",
  "KFAS stoch", "MARSS em stoch", "MARSS bfgs stoch"
)
colnames(vals) <- c("Q", "R1", "Rcov", "R2", "logLik")
round(vals, digits = 5)


###################################################
### code chunk number 63: Cs705_globaltemp
###################################################
out_temp <- KFS(kfas_temp_stoch$model)
df <- data.frame(
  t = as.vector(time(coef(out_temp))),
  KFAS = coef(out_temp),
  `MARSS BFGS` = tsSmooth(marss_temp_stoch_bfgs, type = "xtT")$.estimate,
  `MARSS EM` = tsSmooth(marss_temp_stoch_bfgs, type = "xtT")$.estimate
)
df <- pivot_longer(df, c(KFAS, MARSS.BFGS, MARSS.EM), names_to = "model", values_to = "value")
ggplot(df, aes(x = t, y = value)) +
  geom_line() +
  facet_wrap(~model)


###################################################
### code chunk number 64: Cs706_globaltemp
###################################################
mod.list <- list(
  Z = matrix(kfas_temp_stoch$model$Z, ncol = 1),
  R = kfas_temp_stoch$model$H[, , 1],
  U = matrix(0),
  A = matrix(0, 2, 1),
  Q = matrix(kfas_temp_stoch$model$Q[, , 1]),
  x0 = kfas_temp_stoch$model$a1,
  V0 = kfas_temp_stoch$model$P1,
  tinitx = 1
)
marss_test <- MARSS(t(GlobalTemp), model = mod.list)


###################################################
### code chunk number 65: Cs707_globaltemp
###################################################
kfs <- KFS(kfas_temp_stoch$model, smoothing = "disturbance")
resid_kfas <- rstandard(kfs,
  type = "state",
  standardization_type = "cholesky"
)
resid_marss <- residuals(marss_temp_stoch_bfgs,
  type = "tT",
  standardization = "Block.Cholesky"
)
resid_marss <- subset(resid_marss, name == "state")
resid_test <- residuals(marss_test,
  type = "tT",
  standardization = "Block.Cholesky"
)
resid_test <- subset(resid_test, name == "state")

df <- data.frame(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas),
  MARSS.test = resid_test$.std.resids
)
df$diff.est <- df$MARSS - df$KFAS
df$diff.id <- df$MARSS.test - df$KFAS
df$t <- as.vector(time(kfas_temp_stoch$model$y))
df$name <- "state"
df1 <- pivot_longer(df, c(MARSS, MARSS.test, KFAS, diff.est, diff.id), names_to = "model", values_to = "value")

kfs <- KFS(kfas_temp_stoch$model)
resid_kfas <- residuals(kfs, type = "pearson")
resid_marss <- MARSSresiduals(marss_temp_stoch_bfgs, type = "tT")
resid_test <- MARSSresiduals(marss_test, type = "tT")
df <- rbind(
  cbind(as.data.frame(t(resid_marss$model.residuals)), model = "MARSS") %>% pivot_longer(c("HL", "Folland")),
  cbind(as.data.frame(t(resid_test$model.residuals)), model = "MARSS.test") %>% pivot_longer(c("HL", "Folland")),
  cbind(as.data.frame(resid_kfas), model = "KFAS") %>% pivot_longer(c("HL", "Folland")),
  cbind(as.data.frame(t(resid_marss$model.residuals) - resid_kfas), model = "diff.est") %>% pivot_longer(c("HL", "Folland")),
  cbind(as.data.frame(t(resid_test$model.residuals) - resid_kfas), model = "diff.id") %>% pivot_longer(c("HL", "Folland"))
)
df$t <- rep(as.vector(time(kfas_temp_stoch$model$y)), 2 * 5)

df <- rbind(df, df1[, colnames(df)])

ggplot(subset(df, model %in% c("diff.est", "diff.id")), aes(x = t, y = value, col = model)) +
  geom_line(na.rm = TRUE) +
  facet_wrap(~name) +
  xlab("") +
  ggtitle("Difference in residuals KFAS vs MARSS")


###################################################
### code chunk number 66: Cs708_globaltemp (eval = FALSE)
###################################################
## # test; these should be identical
## kfas_test <- kfas_temp_stoch
## mod.list <- list(
##   Z = matrix(1, 2, 1),
##   R = kfas_test$model$H[, , 1],
##   U = matrix(0),
##   A = matrix(0, 2, 1),
##   Q = matrix(kfas_test$model$Q[, , 1]),
##   tinitx = 1
## )
## mod.list$x0 <- matrix(0)
## mod.list$V0 <- kfas_test$model$P1
## marss_test <- MARSS(t(GlobalTemp), model = mod.list)
## kfas_test_mod <- MARSSkfas(marss_test,
##   return.kfas.model = TRUE,
##   return.lag.one = FALSE
## )$kfas.model
## 
## kfs <- KFS(kfas_test_mod, smoothing = "disturbance")
## test <- cbind(b = kfas_test_mod$Q[1, 1, 1] - kfs$V_eta[1, 1, ], a = MARSSresiduals(marss_test, type = "tT")$var.residuals[3, 3, ])
## test <- as.data.frame(test)
## test$diff <- test$b - test$a
## head(test)
## tail(test)
## 
## MARSSkfas(marss_test)$VtT[, , 1]
## MARSSkfss(marss_test)$VtT[, , 1]
## 
## 
## var.EytT_fit <-
##   fitted(marss_test, type = "ytT", interval = "confidence")$.se^2
## cbind(V_mu = KFS(kfas_test$model)$V_mu[1, 1, ], fitted = var.EytT_fit)
## 
## kfs <- KFS(kfas_test_mod, smoothing = "disturbance")
## resid_kfas <- rstandard(kfs,
##   type = "state",
##   standardization_type = "marginal"
## )
## resid_marss <- residuals(marss_test,
##   type = "tT",
##   standardization = "marginal"
## )
## resid_marss <- subset(resid_marss, name == "state")
## df <- data.frame(
##   MARSS = resid_marss$.std.resids,
##   KFAS = as.vector(resid_kfas),
##   diff = resid_marss$.std.resids - as.vector(resid_kfas)
## )
## head(df)
## 
## resid_kfas <- residuals(kfs, type = "state")
## resid_marss <- residuals(marss_test, type = "tT")
## resid_marss <- subset(resid_marss, name == "state")
## df <- cbind(
##   MARSS = resid_marss$.resids,
##   KFAS = as.vector(resid_kfas),
##   diff = resid_marss$.resids - as.vector(resid_kfas)
## )
## head(df)
## 
## 
## kfs <- KFS(kfas_temp_stoch$model)
## resid_kfas <- residuals(kfs, type = "pearson")
## resid_marss <- MARSSresiduals(marss_test, type = "tT")
## df <- cbind(
##   as.data.frame(t(resid_marss$model.residuals)),
##   as.data.frame(resid_kfas),
##   as.data.frame(t(resid_marss$model.residuals) - resid_kfas)
## )


