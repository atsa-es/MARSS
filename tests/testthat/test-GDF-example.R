skip_on_cran()

library(MARSS)
library(lubridate)
library(tidyverse)
library("tseries")
library("dplyr")
library(quantmod)

# 2) Get data from Quantmod and prepare data frame for estimation

# getSymbols('GDPC1',src='FRED')
# getSymbols('PAYEMS',from = "1947-01-01",src='FRED')
load("GDP.Rdata")

GDP <- data.frame(date = index(GDPC1), coredata(GDPC1))
Emp <- data.frame(date = index(PAYEMS), coredata(PAYEMS))
Emp <- Emp %>% filter(date >= as.Date("1947-01-01") & date <= as.Date("2020-06-01"))

Emp$PAYEMS <- as.numeric(Emp$PAYEMS)
Emp <- Emp %>% mutate(rate = PAYEMS / lag(PAYEMS, 1) - 1)
Emp <- Emp %>% mutate(norm_rate = scale(rate, center = TRUE, scale = TRUE))
Emp <- select(Emp, -c(rate))

GDP <- GDP %>% mutate(rate = GDPC1 / lag(GDPC1, 1) - 1)
GDP <- GDP %>% mutate(norm_rate = scale(rate, center = TRUE, scale = TRUE))
GDP <- select(GDP, -c(rate, GDPC1))

months <- lapply(X = GDP$date, FUN = seq.Date, by = "month", length.out = 3)
months <- data.frame(date = do.call(what = c, months))

m_GDP <- left_join(x = months, y = GDP, by = "date")

df <- cbind(m_GDP, Emp$norm_rate)

names(df) <- c("date", "S01_GDP", "S02_Emp")

df_marss <- df %>% gather(key = "serie", value = "value", -date)

df_marss <- df_marss %>% spread(key = date, value = value)

df_marss$serie <- NULL

df_marss <- as.matrix(df_marss)

# 3) Define State Space matrices

# Matrix Z

Z <- matrix(list(
  "0.33*z1", "z2",
  "0.67*z1", 0,
  "z1", 0,
  "0.67*z1", 0,
  "0.33*z1", 0,
  1 / 3, 0,
  2 / 3, 0,
  1, 0,
  2 / 3, 0,
  1 / 3, 0,
  0, 1,
  0, 0
), 2, 12)

m <- nrow(Z)
p <- ncol(Z)

# Matrix R

R <- matrix(list(0), m, m)

# Matrix B

B <- matrix(list(
  "b1", 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  "b2", 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "b6", 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "b7", 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b11", 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b12", 0
), 12, 12)


# Matrix Q

Q <- matrix(list(
  "q1", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "q6", 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "q11", 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
), 12, 12)

x0 <- matrix(0, p, 1)
A <- matrix(0, (length(df) - 1), 1)
U <- matrix(0, p, 1)
V0 <- 5 * diag(1, p)
U <- matrix(0, p, 1)

# 4) Estimation

# Define model

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 0)

# Estimation

kf_ss <- try(MARSS(df_marss, model = model.gen, method = "BFGS", silent = TRUE), silent = TRUE)

test_that("GDF example for numerical stabilty", {
  expect_true(!inherits(kf_ss, "try-error"))
})

test_that("GDF example for numerical stabilty", {
  expect_true(kf_ss$convergence == 54)
})

test_that("GDF example for numerical stabilty", {
  expect_true(all.equal(kf_ss$logLik, -86.3201321))
})

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 1)
kf_ss <- try(MARSS(df_marss, model = model.gen, method = "BFGS", silent = TRUE), silent = TRUE)

test_that("GDF example for numerical stabilty", {
  expect_true(!inherits(kf_ss, "try-error"))
})

test_that("GDF example for numerical stabilty", {
  expect_true(kf_ss$convergence == 54)
})

test_that("GDF example for numerical stabilty", {
  expect_true(all.equal(kf_ss$logLik, -991.2024293))
})

# THIS One should work; scaling is different

GDP <- data.frame(date = index(GDPC1), coredata(GDPC1))
Emp <- data.frame(date = index(PAYEMS), coredata(PAYEMS))
Emp <- Emp %>% filter(date >= as.Date("1947-01-01") & date <= as.Date("2020-06-01"))

Emp$PAYEMS <- as.numeric(Emp$PAYEMS)
Emp <- Emp %>% mutate(rate = PAYEMS / lag(PAYEMS, 1) - 1)
GDP <- GDP %>% mutate(rate = GDPC1 / lag(GDPC1, 1) - 1)
GDP <- select(GDP, -c(GDPC1))

months <- lapply(X = GDP$date, FUN = seq.Date, by = "month", length.out = 3)
months <- data.frame(date = do.call(what = c, months))

m_GDP <- left_join(x = months, y = GDP, by = "date")

df <- cbind(m_GDP, Emp$rate)

names(df) <- c("date", "S01_GDP", "S02_Emp")

df_marss <- df %>% gather(key = "serie", value = "value", -date)

df_marss <- df_marss %>% spread(key = date, value = value)

df_marss$serie <- NULL

df_marss <- as.matrix(df_marss)


# Define model

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 0)
kf_ss <- try(MARSS(df_marss, model = model.gen, method = "BFGS", silent = TRUE), silent = TRUE)

# Estimation

test_that("GDF works example for numerical stabilty", {
  expect_true(!inherits(kf_ss, "try-error"))
})

test_that("GDF works example for numerical stabilty", {
  expect_true(kf_ss$convergence == 0)
})

test_that("GDF works example for numerical stabilty", {
  expect_true(all.equal(kf_ss$logLik, 4260.399772))
})

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 1)
kf_ss <- try(MARSS(df_marss, model = model.gen, method = "BFGS", silent = TRUE), silent = TRUE)

test_that("GDF works example for numerical stabilty", {
  expect_true(!inherits(kf_ss, "try-error"))
})

test_that("GDF works example for numerical stabilty", {
  expect_true(kf_ss$convergence == 0)
})

test_that("GDF works example for numerical stabilty", {
  expect_true(all.equal(kf_ss$logLik, 4301.624359))
})

### More difficult

# getSymbols('GDPC1',src='FRED')
# getSymbols('PAYEMS',from = "1947-01-01",src='FRED')
# getSymbols('INDPRO',from = "1947-01-01",src='FRED')
# getSymbols('RPI',from = "1947-01-01",src='FRED')
# save(GDPC1, PAYEMS, INDPRO, RPI, file="tests/testthat/GDP.RData")
load("GDP.Rdata")

GDP <- data.frame(date = index(GDPC1), coredata(GDPC1))
Emp <- data.frame(date = index(PAYEMS), coredata(PAYEMS))
Indpr <- data.frame(date = index(INDPRO), coredata(INDPRO))
Inc <- data.frame(date = index(RPI), coredata(RPI))

Emp <- Emp %>% filter(date >= as.Date("1947-01-01") & date <= as.Date("2020-06-01"))
Indpr <- Indpr %>% filter(date >= as.Date("1947-01-01") & date <= as.Date("2020-06-01"))
Inc <- Inc %>% filter(date >= as.Date("1947-01-01") & date <= as.Date("2020-06-01"))

names(Inc) <- c("date", "Inc")

Inc_aux <- data.frame(seq.Date(as.Date("1947-01-01"), as.Date("1958-12-01"), by = "month"))
Inc_aux$RPI <- NA

names(Inc_aux) <- c("date", "Inc")

Inc <- rbind(Inc_aux, Inc)

Emp$PAYEMS <- as.numeric(Emp$PAYEMS)
Emp <- Emp %>% mutate(rate = PAYEMS / lag(PAYEMS, 1) - 1)

Indpr$INDPRO <- as.numeric(Indpr$INDPRO)
Indpr <- Indpr %>% mutate(rate = INDPRO / lag(INDPRO, 1) - 1)

Inc$Inc <- as.numeric(Inc$Inc)
Inc <- Inc %>% mutate(rate = Inc / lag(Inc, 1) - 1)


GDP <- GDP %>% mutate(rate = GDPC1 / lag(GDPC1, 1) - 1)
GDP <- select(GDP, -c(GDPC1))

months <- lapply(X = GDP$date, FUN = seq.Date, by = "month", length.out = 3)
months <- data.frame(date = do.call(what = c, months))

m_GDP <- left_join(x = months, y = GDP, by = "date")

# Data frame for estimation

df <- cbind(m_GDP, Emp$rate, Indpr$rate, Inc$rate)

names(df) <- c("date", "S01_GDP", "S02_Emp", "S03_Indpr", "S04_Inc")

# 3) Model 1: one quarterly series and two monthly series

df_marss <- select(df, -c(S04_Inc))

df_marss <- df_marss %>% gather(key = "serie", value = "value", -date)

df_marss <- df_marss %>% spread(key = date, value = value)

df_marss$serie <- NULL

df_marss <- as.matrix(df_marss)

df_marss <- zscore(df_marss, mean.only = TRUE)


# Matrix Z
#
Z <- matrix(list(
  "0.33*z1", "z2", "z3",
  "0.67*z1", 0, 0,
  "z1", 0, 0,
  "0.67*z1", 0, 0,
  "0.33*z1", 0, 0,
  1 / 3, 0, 0,
  2 / 3, 0, 0,
  1, 0, 0,
  2 / 3, 0, 0,
  1 / 3, 0, 0,
  0, 1, 0,
  0, 0, 0,
  0, 0, 1,
  0, 0, 0
), 3, 14)

m <- nrow(Z)
p <- ncol(Z)

# Matrix R

R <- matrix(list(0), m, m)

# Matrix B

B <- matrix(list(
  "b1", 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  "b2", 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "b6", 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "b7", 0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b11", 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b12", 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b13", 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b14", 0
), 14, 14)


# Matrix Q

Q <- matrix(list(
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "q6", 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "q11", 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "q13", 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
), 14, 14)

Q[1, 1] <- "q1"

# Rest of matrices

x0 <- matrix(0, p, 1)
A <- matrix(0, (length(df) - 2), 1)
U <- matrix(0, p, 1)
V0 <- 5 * diag(1, p)
U <- matrix(0, p, 1)


# Estimation

# Define model

model.gen_1 <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 1)

# Estimation

kf_ss <- try(MARSS(df_marss, model = model.gen_1, method = "BFGS", silent = TRUE), silent = TRUE)

test_that("GDF example 3 for numerical stabilty", {
  expect_true(!inherits(kf_ss, "try-error"))
})

test_that("GDF example 3 for numerical stabilty", {
  expect_true(kf_ss$convergence == 0)
})

test_that("GDF example 3 for numerical stabilty", {
  expect_true(all.equal(kf_ss$logLik, 7380.469127))
})

# 4) Model 2: one quarterly series and three monthly series

df_marss <- df %>% gather(key = "serie", value = "value", -date)

df_marss <- df_marss %>% spread(key = date, value = value)

df_marss$serie <- NULL

df_marss <- as.matrix(df_marss)

df_marss <- zscore(df_marss, mean.only = TRUE)

# Matrix Z
#
Z <- matrix(list(
  "0.33*z1", "z2", "z3", "z4",
  "0.67*z1", 0, 0, 0,
  "z1", 0, 0, 0,
  "0.67*z1", 0, 0, 0,
  "0.33*z1", 0, 0, 0,
  1 / 3, 0, 0, 0,
  2 / 3, 0, 0, 0,
  1, 0, 0, 0,
  2 / 3, 0, 0, 0,
  1 / 3, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 0,
  0, 0, 0, 1,
  0, 0, 0, 0
), 4, 16)

m <- nrow(Z)
p <- ncol(Z)

# Matrix R

R <- matrix(list(0), m, m)

# Matrix B

B <- matrix(list(
  "b1", 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  "b2", 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "b6", 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "b7", 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b11", 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b12", 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b13", 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b14", 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b15", 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "b16", 0
), 16, 16)


# Matrix Q

Q <- matrix(list(
  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, "q6", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "q11", 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "q13", 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "q15", 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
), 16, 16)

Q[1, 1] <- "q1"

# Rest of matrices

x0 <- matrix(0, p, 1)
A <- matrix(0, (length(df) - 1), 1)
U <- matrix(0, p, 1)
V0 <- 5 * diag(1, p)
U <- matrix(0, p, 1)

# Estimation

# Define model

model.gen_2 <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 1)

# Estimation

kf_ss <- try(MARSS(df_marss, model = model.gen_2, method = "BFGS", silent = TRUE), silent = TRUE)

test_that("GDF example 4 for numerical stabilty", {
  expect_true(!inherits(kf_ss, "try-error"))
})

test_that("GDF example 4 for numerical stabilty", {
  expect_true(kf_ss$convergence == 0)
})

test_that("GDF example 4 for numerical stabilty", {
  expect_true(all.equal(kf_ss$logLik, 10014.90035))
})
