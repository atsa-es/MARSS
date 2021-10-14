model.list <- list()

name <- "RW_1D_a"
set.seed(123)
dat <- cumsum(rnorm(20))
mod.list <- list(tinitx = 1, U = "zero", R = "unconstrained", Q = "unconstrained", B = "unconstrained", x0 = matrix(dat[1] * 1.0001))
con.list <- NULL
m <- 1
n <- 1
model.list[[name]] <- list(data = dat, model = mod.list, control = con.list, name = name, m = m, n = n, kfss = TRUE)

name <- "RW_1D_b"
mod.list <- list(tinitx = 0, U = "zero", R = "unconstrained", Q = "unconstrained", B = "unconstrained")
model.list[[name]] <- list(data = dat, model = mod.list, control = con.list, name = name, m = m, n = n, kfss = TRUE)

# Little harder model
namebase <- "HarborSealWA234"
dat <- t(harborSealWA)
dat <- dat[2:4, ] # remove the year row
con.list <- NULL
m <- 3
n <- 3

for (Q in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")) {
  for (B in c("identity", "diagonal and unequal")) {
    for (R in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")) {
      if (B == "diagonal and unequal" && (is.list(Q) || is.list(R))) next
      mod.list <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = dat[, 1, drop = FALSE] * 1.1)
      if (B != "identity" && R != "zero") mod.list$tinitx <- 1
      if (Q == "zero" && R == "zero") next
      if (Q == "zero" && B != "identity") next
      name <- paste(namebase, Q, R, B, sep = "-")
      model.list[[name]] <- list(data = dat, model = mod.list, control = con.list, name = name, m = m, n = n, kfss = TRUE)
    }
  }
}

# test B unconstrained

Q <- "diagonal and unequal"
B <- "unconstrained"
R <- "diagonal and unequal"
name <- paste(namebase, Q, R, B, sep = "-")
mod.list <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = "unequal", tinitx = 1)
model.list[[name]] <- list(data = dat, model = mod.list, control = con.list, name = name, m = m, n = n, kfss = TRUE)

# test Q with some zeros

Q <- ldiag(list("q1", 0, "q2"))
B <- "identity"
R <- "diagonal and equal"
name <- paste(namebase, "QwZero", R, B, sep = "-")
mod.list <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = "unequal", tinitx = 1)
model.list[[name]] <- list(data = dat, model = mod.list, control = con.list, name = name, m = m, n = n, kfss = TRUE)

# test R with some zeros

R <- ldiag(list("r1", 0, "r2"))
B <- "identity"
Q <- "diagonal and equal"
name <- paste(namebase, Q, "RwZeros", B, sep = "-")
mod.list <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = "unequal", tinitx = 1)
model.list[[name]] <- list(data = dat, model = mod.list, control = con.list, name = name, m = m, n = n, kfss = TRUE)

# Wonky model; this is simple version of the GDP test
# 1) Define some data

df_marss <- matrix(NA, 2, 10)
df_marss[1, ] <- c(NA, NA, NA, -0.002666915, NA, NA, -0.002064963, NA, NA, 0.01564208)
df_marss[2, ] <- c(NA, 0.0005053405, 0.001147921, -0.002476667, 0.003195476, 0.003941519, -0.001529331, 0.004960794, 0.005527753, 0.004705563)

# 2) Define State Space matrices

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

# Rest of matrices

x0 <- matrix(0, p, 1)
A <- matrix(0, m, 1)
U <- matrix(0, p, 1)
V0 <- 5 * diag(1, p)
U <- matrix(0, p, 1)

# Define model
model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 0)
con.list <- NULL
name <- "GDP1"
model.list[[name]] <- list(data = df_marss, model = model.gen, control = con.list, name = name, m = p, n = m, kfss = FALSE)

# More GDP Models

library(lubridate)
library(tidyverse)
library("tseries")
library("dplyr")
library(quantmod)

# 2) Get data from Quantmod and prepare data frame for estimation

getSymbols("GDPC1", src = "FRED")
getSymbols("PAYEMS", from = "1947-01-01", src = "FRED")

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
m_GDP <- m_GDP %>% filter(date >= as.Date("1947-01-01") & date <= as.Date("2020-06-01"))

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


# This one throws an error

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 0)
con.list <- NULL
name <- "GDP3-tinit0"
model.list[[name]] <- list(data = df_marss, model = model.gen, control = con.list, name = name, m = p, n = m, kfss = FALSE)

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 1)
name <- "GDP3-tinit1"
model.list[[name]] <- list(data = df_marss, model = model.gen, control = con.list, name = name, m = p, n = m, kfss = FALSE)

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
m_GDP <- m_GDP %>% filter(date >= as.Date("1947-01-01") & date <= as.Date("2020-06-01"))

df <- cbind(m_GDP, Emp$rate)
names(df) <- c("date", "S01_GDP", "S02_Emp")
df_marss <- df %>% gather(key = "serie", value = "value", -date)
df_marss <- df_marss %>% spread(key = date, value = value)
df_marss$serie <- NULL
df_marss <- as.matrix(df_marss)

# Define model

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 0)
name <- "GDP2-tinit0"
model.list[[name]] <- list(data = df_marss, model = model.gen, control = con.list, name = name, m = p, n = m, kfss = FALSE)

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 1)
name <- "GDP2-tinit1"
model.list[[name]] <- list(data = df_marss, model = model.gen, control = con.list, name = name, m = p, n = m, kfss = FALSE)


### More difficult

getSymbols("GDPC1", src = "FRED")
getSymbols("PAYEMS", from = "1947-01-01", src = "FRED")
getSymbols("INDPRO", from = "1947-01-01", src = "FRED")
getSymbols("RPI", from = "1947-01-01", src = "FRED")

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
m_GDP <- m_GDP %>% filter(date >= as.Date("1947-01-01") & date <= as.Date("2020-06-01"))

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


# Define model
model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 1)
name <- "GDP4-tinit1"
model.list[[name]] <- list(data = df_marss, model = model.gen, control = con.list, name = name, m = p, n = (length(df) - 2), kfss = FALSE)


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

model.gen <- list(Z = Z, A = A, R = R, B = B, U = U, Q = Q, x0 = x0, V0 = V0, tinitx = 1)
name <- "GDP5-tinit1"
model.list[[name]] <- list(data = df_marss, model = model.gen, control = con.list, name = name, m = p, n = (length(df) - 1), kfss = FALSE)

## StructTS models

name <- "StructTS-treering-fixed"
y <- window(treering, start = 0, end = 20)
fit1 <- StructTS(y, type = "level")

# Run with fit1 estimates
vy <- var(y, na.rm = TRUE) / 100
mod.list <- list(
  x0 = matrix(y[1]), U = "zero", tinitx = 0,
  Q = matrix(fit1$coef[1]), R = matrix(fit1$coef[2]),
  V0 = matrix(1e+06 * vy)
)
model.list[[name]] <- list(data = as.vector(y), model = mod.list, control = NULL, name = name, m = 1, n = 1, kfss = TRUE)

# Now estimate the parameters
name <- "StructTS-treering"
mod.list <- list(
  x0 = matrix(y[1]), U = "zero", tinitx = 0, V0 = matrix(1e+06 * vy),
  Q = matrix("s2xi"), R = matrix("s2eps")
)
model.list[[name]] <- list(data = as.vector(y), model = mod.list, control = list(allow.degen = FALSE), name = name, m = 1, n = 1, kfss = TRUE)

### Nile with NAs

mod.nile <- list(
  Z = matrix(1), A = matrix(0), R = matrix("r"),
  B = matrix(1), U = matrix(0), Q = matrix("q"),
  tinitx = 1
)

dat <- t(as.matrix(Nile))
rownames(dat) <- "Nile"
name <- "StructTS-Nile-level"
model.list[[name]] <- list(data = dat, model = mod.nile, control = NULL, name = name, m = 1, n = 1, kfss = TRUE)

### Nile with NAs

NileNA <- Nile
NileNA[c(21:40, 61:80)] <- NA
name <- "StructTS-NileNA-level"
model.list[[name]] <- list(data = as.vector(NileNA), model = mod.nile, control = NULL, name = name, m = 1, n = 1, kfss = TRUE)

###################################################
### Global Temp
###################################################

data("GlobalTemp", package = "KFAS")

mod.list <- list(
  Z = matrix(1, 2, 1),
  R = matrix(c("r1", "c", "c", "r2"), 2, 2),
  U = matrix(0),
  A = matrix(0, 2, 1),
  tinitx = 1
)
name <- "StructTS-GlobalTemp"
model.list[[name]] <- list(data = t(GlobalTemp), model = mod.list, control = NULL, name = name, m = 1, n = 2, notes = "BFGS only", kfss = TRUE)

save(model.list, file = "tests/testthat/models.RData")
