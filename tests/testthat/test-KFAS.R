skip_on_cran()

library(MARSS)
library(KFAS)

###################################################
### Nile
###################################################
context("KFAS comparison Nile")

model_Nile <- SSModel(Nile ~ SSMtrend(
  degree = 1,
  Q = list(matrix(NA))
),
H = matrix(NA)
)
kinits <- c(log(var(Nile)), log(var(Nile)))
model_Nile_stoch <- model_Nile
model_Nile_stoch$a1[1, 1] <- 0
model_Nile_stoch$P1[1, 1] <- 1000
model_Nile_stoch$P1inf[1, 1] <- 0
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_kfas_stoch <- fitSSM(model_Nile_stoch, kinits, method = "BFGS")
kfs_kfas_stoch <- KFS(fit_kfas_stoch$model)

# MARSS
mod.nile <- list(
  Z = matrix(1), A = matrix(0), R = matrix("r"),
  B = matrix(1), U = matrix(0), Q = matrix("q"),
  tinitx = 1
)

dat <- t(as.matrix(Nile))
rownames(dat) <- "Nile"
inits <- list(Q = matrix(var(Nile)), R = matrix(var(Nile)))
mod.nile.stoch <- mod.nile
mod.nile.stoch$x0 <- fit_kfas_stoch$model$a1
mod.nile.stoch$V0 <- fit_kfas_stoch$model$P1
fit_em_stoch <- MARSS(dat, model = mod.nile.stoch, silent = TRUE)
fit_bfgs_stoch <- MARSS(dat,
  model = mod.nile.stoch, inits = inits,
  method = "BFGS", silent = TRUE
)

marss_kfas_model <- MARSSkfas(fit_em_stoch,
  return.kfas.model = TRUE,
  return.lag.one = FALSE
)$kfas.model
marss_kfas_model$Q[1, 1, 1] <- NA
marss_kfas_model$H[1, 1, 1] <- NA
kinits <- c(log(var(Nile)), log(var(Nile)))
fit_marss_kfas <- fitSSM(marss_kfas_model, kinits, method = "BFGS")

vals <- rbind(
  c(fit_kfas_stoch$model$Q, fit_kfas_stoch$model$H, -1 * fit_kfas_stoch$optim.out$value),
  c(coef(fit_em_stoch)$Q, coef(fit_em_stoch)$R, logLik(fit_em_stoch)),
  c(coef(fit_bfgs_stoch)$Q, coef(fit_bfgs_stoch)$R, logLik(fit_bfgs_stoch)),
  c(fit_marss_kfas$model$Q, fit_marss_kfas$model$H, -1 * fit_marss_kfas$optim.out$value)
)
rownames(vals) <- c(
  "KFAS",
  "MARSSem",
  "MARSSbfgs",
  "KFASmarss"
)
colnames(vals) <- c("Q", "R", "logLik")
vals <- as.data.frame(vals)

test_that("fits Nile", {
  expect_equivalent(vals["KFAS", ], vals["KFASmarss", ])
})

test_that("fits Nile KFAS vs MARSS", {
  expect_true(all((vals$logLik - max(vals$logLik[1])) / max(vals$logLik[1]) <= 1.293388e-06))
})


###################################################
### Test state filtering
###################################################

fit_kfas <- fit_kfas_stoch
fit_marss <- fit_em_stoch
fit_marss$par$Q[1, 1] <- fit_kfas$model$Q
fit_marss$par$R[1, 1] <- fit_kfas$model$H

kf_kfas <- KFS(fit_kfas$model,
  filtering = "state",
  smoothing = "state", simplify = FALSE
)
kf_marss <- MARSSkfss(fit_marss)

test_that("fits Nile xtt1", {
  expect_equivalent(kf_kfas$a[1:100], kf_marss$xtt1[, 1:100])
})
test_that("fits Nile xtt", {
  expect_equivalent(kf_kfas$att[1:100], kf_marss$xtt[, 1:100])
})
test_that("fits Nile xtT", {
  expect_equivalent(kf_kfas$alphahat[1:100], kf_marss$xtT[, 1:100])
})
test_that("fits Nile Innov", {
  expect_equivalent(kf_kfas$v[1:100], kf_marss$Innov[, 1:100])
})
test_that("fits Nile Sigma", {
  expect_equivalent(kf_kfas$F[1:100], kf_marss$Sigma[1:100])
})
test_that("fits Nile Vtt1", {
  expect_equivalent(kf_kfas$P[1:100], kf_marss$Vtt1[1:100])
})
test_that("fits Nile Vtt", {
  expect_equivalent(kf_kfas$Ptt[1:100], kf_marss$Vtt[1:100])
})

###################################################
### observation filtering tests
###################################################
f_kfas <- KFS(fit_kfas$model,
  filtering = "signal",
  smoothing = "signal", simplify = FALSE
)
f_marss <- MARSSkfss(fit_marss)

n <- 100

ytt1_fit <- fitted(fit_marss, type = "ytt1")$.fitted
ytt1_hatyt <- MARSShatyt(fit_marss, only.kem = FALSE)$ytt1
test_that("fits Nile ytt1", {
  expect_equivalent(f_kfas$m[1:n], ytt1_fit[1:n], ytt1_hatyt[1:n])
})

var.Eytt1_fit <-
  fitted(fit_marss, type = "ytt1", interval = "confidence")$.se^2
var.Eytt1_hatyt <-
  MARSShatyt(fit_marss, only.kem = FALSE)$var.Eytt1
test_that("fits Nile Eytt1", {
  expect_equivalent(f_kfas$P_mu[1:n], var.Eytt1_fit[1:n])
})

ytT_fit <- fitted(fit_marss, type = "ytT")$.fitted
ytT_hatyt <- MARSShatyt(fit_marss)$ytT
test_that("fits Nile ytT", {
  expect_equivalent(f_kfas$muhat[1:n], ytT_fit[1:n])
})

var.EytT_fit <-
  fitted(fit_marss, type = "ytT", interval = "confidence")$.se^2
var.EytT_hatyt <-
  MARSShatyt(fit_marss, only.kem = FALSE)$var.EytT
test_that("fits Nile var.EytT", {
  expect_equivalent(f_kfas$V_mu[1:n], var.EytT_fit[1:n])
})

###################################################
### Predicting
###################################################

conf_kfas <- predict(fit_kfas$model,
  interval = "confidence",
  se.fit = TRUE
)
conf_marss1 <- fitted(fit_marss, type = "ytT", interval = "confidence")

conf_marss2 <- predict(fit_marss,
  type = "ytT",
  interval = "confidence", level = 0.95
)
class(conf_kfas) <- "matrix"
test_that("fits Nile CI", {
  expect_equivalent(conf_kfas, as.matrix(conf_marss1[, c(".fitted", ".conf.low", ".conf.up", ".se")]), as.matrix(conf_marss2$pred[, c("estimate", "Lo 95", "Hi 95", "se")]))
})

pred_kfas <- predict(fit_kfas$model,
  interval = "prediction", se.fit = TRUE
)
pred_marss1 <- fitted(fit_marss, type = "ytT", interval = "prediction")
pred_marss2 <- predict(fit_marss,
  type = "ytT",
  interval = "prediction", level = 0.95
)
class(pred_kfas) <- "matrix"
test_that("fits Nile PI", {
  expect_equivalent(pred_kfas, as.matrix(pred_marss1[, c(".fitted", ".lwr", ".upr", ".sd")]), as.matrix(pred_marss2$pred[, c("estimate", "Lo 95", "Hi 95", "se")]))
})

conf_kfas_t1 <- predict(fit_kfas$model,
  interval = "confidence",
  se.fit = TRUE, filtered = TRUE
)
class(conf_kfas_t1) <- "matrix"
conf_marss1_t1 <- fitted(fit_marss, type = "ytt1", interval = "confidence")
conf_marss2_t1 <- predict(fit_marss,
  type = "ytt1",
  interval = "confidence", level = 0.95
)
test_that("fits Nile CI t1", {
  expect_equivalent(conf_kfas_t1, as.matrix(conf_marss1_t1[, c(".fitted", ".conf.low", ".conf.up", ".se")]), as.matrix(conf_marss2_t1$pred[, c("estimate", "Lo 95", "Hi 95", "se")]))
})

pred_kfas_t1 <- predict(fit_kfas$model,
  interval = "prediction",
  se.fit = TRUE, filtered = TRUE
)

class(pred_kfas_t1) <- "matrix"
pred_marss1_t1 <- fitted(fit_marss, type = "ytt1", interval = "prediction")
pred_marss2_t1 <- predict(fit_marss,
  type = "ytt1",
  interval = "prediction", level = 0.95
)
test_that("fits Nile PI t1", {
  expect_equivalent(pred_kfas_t1, as.matrix(pred_marss1_t1[, c(".fitted", ".lwr", ".upr", ".sd")]), as.matrix(pred_marss2_t1$pred[, c("estimate", "Lo 95", "Hi 95", "se")]))
})

###################################################
### Residuals
###################################################

kfs <- KFS(fit_kfas$model)
resid_kfas <- residuals(kfs, type = "recursive")
resid_marss <- residuals(fit_marss,
  type = "tt1",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "model")
test_that("fits Nile model residuals recursive", {
  expect_equivalent(resid_marss$.resids, as.vector(resid_kfas))
})

resid_kfas <- rstandard(kfs,
  type = "recursive",
  standardization_type = "marginal"
)
resid_marss <- residuals(fit_marss,
  type = "tt1",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "model")
test_that("fits Nile model residuals marginal", {
  expect_equivalent(resid_marss$.std.resids[1:(n - 1)], as.vector(resid_kfas)[1:(n - 1)])
})

resid_kfas <- rstandard(kfs,
  type = "recursive",
  standardization_type = "cholesky"
)
resid_marss <- residuals(fit_marss,
  type = "tt1",
  standardization = "Cholesky"
)
resid_marss <- subset(resid_marss, name == "model")
test_that("fits Nile model residuals cholesky", {
  expect_equivalent(resid_marss$.std.resids[1:(n - 1)], as.vector(resid_kfas)[1:(n - 1)])
})

resid_kfas <- residuals(kfs, type = "pearson")
resid_marss <- residuals(fit_marss, type = "tT")
resid_marss <- subset(resid_marss, name == "model")
test_that("fits Nile model residuals pearson", {
  expect_equivalent(resid_marss$.resids, as.vector(resid_kfas))
})

resid_kfas <- rstandard(kfs,
  type = "pearson",
  standardization_type = "marginal"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "model")
test_that("fits Nile model residuals pearson marginal", {
  expect_equivalent(resid_marss$.std.resids[1:(n - 1)], as.vector(resid_kfas)[1:(n - 1)])
})

resid_kfas <- rstandard(kfs,
  type = "pearson",
  standardization_type = "cholesky"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "Cholesky"
)
resid_marss <- subset(resid_marss, name == "model")
test_that("fits Nile model residuals pearson cholesky", {
  expect_equivalent(resid_marss$.std.resids[1:(n - 1)], as.vector(resid_kfas)[1:(n - 1)])
})

resid_kfas <- residuals(kfs, type = "response")
resid_marss <- residuals(fit_marss, type = "tT")
resid_marss <- subset(resid_marss, name == "model")
test_that("fits Nile model residuals response", {
  expect_equivalent(resid_marss$.resids, as.vector(resid_kfas))
})

kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
resid_kfas <- residuals(kfs, type = "state")
resid_marss <- residuals(fit_marss, type = "tT")
resid_marss <- subset(resid_marss, name == "state")
test_that("fits Nile state residuals", {
  expect_equivalent(resid_marss$.resids[1:(n - 1)], as.vector(resid_kfas)[1:(n - 1)])
})

resid_kfas <- rstandard(kfs,
  type = "state",
  standardization_type = "marginal"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "marginal"
)
resid_marss <- subset(resid_marss, name == "state")
test_that("fits Nile state residuals marginal", {
  expect_equivalent(resid_marss$.std.resids[1:(n - 1)], as.vector(resid_kfas)[1:(n - 1)])
})

resid_kfas <- rstandard(kfs,
  type = "state",
  standardization_type = "cholesky"
)
resid_marss <- residuals(fit_marss,
  type = "tT",
  standardization = "Block.Cholesky"
)
resid_marss <- subset(resid_marss, name == "state")
test_that("fits Nile state residuals block cholesky", {
  expect_equivalent(resid_marss$.std.resids[1:(n - 1)], as.vector(resid_kfas)[1:(n - 1)])
})

kfs <- KFS(fit_kfas$model, smoothing = "disturbance")
test <- cbind(
  b = fit_kfas$model$Q[1, 1, 1] - kfs$V_eta[1, 1, ],
  a = MARSSresiduals(fit_marss, type = "tT")$var.residuals[2, 2, ]
)
test <- as.data.frame(test)
test_that("fits Nile variance of residuals", {
  expect_equivalent(test$b[1:(n - 1)], test$a[1:(n - 1)])
})

###################################################
### Nile with NAs
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
vals <- rbind(
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
vals <- as.data.frame(vals)

test_that("fits NileNA KFAS vs MARSS", {
  expect_true(all(abs((vals$logLik - max(vals$logLik[1])) / max(vals$logLik[1])) <= 1.5e-06))
})

conf_kfas_NA <-
  predict(fit_kfas_NA$model, interval = "confidence", filtered = FALSE)
conf_marss_NA <-
  predict(fit_marss_NA, interval = "confidence", type = "ytT", level = 0.95)$pred
df1 <- as.data.frame(conf_kfas_NA)
df2 <- conf_marss_NA[, c("estimate", "Lo 95", "Hi 95")]
colnames(df2) <- colnames(df1)
test_that("fits NileNA KFAS vs MARSS CI", {
  expect_true(all(abs((df1 - df2) / df1) < 0.012))
})

df1 <- data.frame(
  smooth = as.vector(fitted(fit_kfas_NA$model)),
  one.step.ahead = as.vector(fitted(fit_kfas_NA$model, filtered = TRUE))
)
df2 <- data.frame(
  smooth = fitted(fit_marss_NA, type = "ytT")$.fitted,
  one.step.ahead = fitted(fit_marss_NA, type = "ytt1")$.fitted
)
test_that("fits NileNA KFAS vs MARSS fitted", {
  expect_true(all(na.omit(abs((df1 - df2) / df1)) < 0.009))
})

###################################################
### Global Temp
###################################################
context("KFAS GlobalTemp fits")

data("GlobalTemp")
model_temp <- SSModel(GlobalTemp ~ SSMtrend(1, Q = NA, type = "common"),
  H = matrix(NA, 2, 2)
)
kinits <- chol(cov(GlobalTemp))[c(1, 4, 3)]
kinits <- c(0.5 * log(0.1), log(kinits[1:2]), kinits[3])
model_temp_stoch <- model_temp
model_temp_stoch$a1[1, 1] <- 0
model_temp_stoch$P1[1, 1] <- 1000 * max(diag(var(GlobalTemp)))
model_temp_stoch$P1inf[1, 1] <- 0
kfas_temp_stoch <- fitSSM(model_temp_stoch, kinits, method = "BFGS")

mod.list <- list(
  Z = matrix(1, 2, 1),
  R = matrix(c("r1", "c", "c", "r2"), 2, 2),
  U = matrix(0),
  A = matrix(0, 2, 1),
  tinitx = 1
)
mod.list$x0 <- kfas_temp_stoch$model$a1
mod.list$V0 <- kfas_temp_stoch$model$P1
marss_temp_stoch_em <- MARSS(t(GlobalTemp), model = mod.list, silent = TRUE)
# use inits from a short run of EM algorithm
inits <- MARSS(t(GlobalTemp),
  model = mod.list, control = list(maxit = 20),
  silent = TRUE
)
marss_temp_stoch_bfgs <- MARSS(t(GlobalTemp), model = mod.list, inits = inits, method = "BFGS", silent = TRUE)

vals <- rbind(
  c(kfas_temp_stoch$model$Q, kfas_temp_stoch$model$H[c(1, 2, 4)], -1 * kfas_temp_stoch$optim.out$value),
  c(coef(marss_temp_stoch_em)$Q, coef(marss_temp_stoch_em)$R, logLik(marss_temp_stoch_em)),
  c(coef(marss_temp_stoch_bfgs)$Q, coef(marss_temp_stoch_bfgs)$R, logLik(marss_temp_stoch_bfgs))
)
rownames(vals) <- c("KFAS stoch", "MARSS em stoch", "MARSS bfgs stoch")
colnames(vals) <- c("Q", "R1", "Rcov", "R2", "logLik")
vals <- as.data.frame(vals)
test_that("fits NileNA KFAS vs MARSS fitted", {
  expect_true((vals$logLik[1] - vals$logLik[3]) / vals$logLik[1] < 1.25e-07)
})
test_that("fits NileNA KFAS vs MARSS fitted", {
  expect_true((vals$logLik[1] - vals$logLik[2]) / vals$logLik[1] < 9.3e-05)
})

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
marss_test <- MARSS(t(GlobalTemp), model = mod.list, silent = TRUE)

out_temp <- KFS(kfas_temp_stoch$model)
df <- data.frame(
  t = as.vector(time(coef(out_temp))),
  KFAS = coef(out_temp),
  `MARSS BFGS` = tsSmooth(marss_temp_stoch_bfgs, type = "xtT")$.estimate,
  `MARSS EM` = tsSmooth(marss_temp_stoch_em, type = "xtT")$.estimate,
  `MARSS test` = tsSmooth(marss_test, type = "xtT")$.estimate
)


test_that("xtT GlobalTemp identical model", {
  expect_equivalent(as.numeric(df$KFAS), df$MARSS.test)
})
test_that("xtT GlobalTemp MARSS em vs KFAS", {
  expect_true(all(abs((df$KFAS - df$MARSS.EM)) < 0.0086))
})
test_that("xtT GlobalTemp MARSS bfgs vs KFAS", {
  expect_true(all(abs((df$KFAS - df$MARSS.BFGS)) < 0.00024))
})

# residuals

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
df <- na.omit(df)

test_that("state block cholesky resids GlobalTemp MARSS equivalent vs KFAS", {
  expect_equal(df$KFAS, df$MARSS.test)
})
test_that("state block cholesky resids GlobalTemp MARSS vs KFAS", {
  expect_true(all(abs(df$KFAS - df$MARSS) < 0.0036))
})

library(tidyr)
kfs <- KFS(kfas_temp_stoch$model)
resid_kfas <- residuals(kfs, type = "pearson")
resid_marss <- MARSSresiduals(marss_temp_stoch_bfgs, type = "tT")
resid_test <- MARSSresiduals(marss_test, type = "tT")
df <- data.frame(
  MARSS = as.data.frame(t(resid_marss$model.residuals)) %>% pivot_longer(c("HL", "Folland")),
  MARSS.test = as.data.frame(t(resid_test$model.residuals)) %>% pivot_longer(c("HL", "Folland")),
  KFAS = as.data.frame(resid_kfas) %>% pivot_longer(c("HL", "Folland"))
)
df <- na.omit(df)
test_that("pearson resids GlobalTemp MARSS equivalent vs KFAS", {
  expect_equal(df$KFAS.value, df$MARSS.test.value)
})
test_that("pearson resids GlobalTemp MARSS vs KFAS", {
  expect_true(all(abs(df$KFAS.value - df$MARSS.value) < 0.000232))
})
