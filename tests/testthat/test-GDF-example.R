skip_on_cran()

library(MARSS)
library(lubridate)
library(tidyverse)
library('tseries')
library('dplyr')
library(quantmod)

# 2) Get data from Quantmod and prepare data frame for estimation

getSymbols('GDPC1',src='FRED')
getSymbols('PAYEMS',from = "1947-01-01",src='FRED')

GDP <-  data.frame(date=index(GDPC1), coredata(GDPC1))
Emp <-  data.frame(date=index(PAYEMS), coredata(PAYEMS))
Emp <- Emp %>%  filter(date>=as.Date("1947-01-01")&date<=as.Date("2020-06-01"))

Emp$PAYEMS <- as.numeric(Emp$PAYEMS)
Emp <- Emp %>% mutate(rate = PAYEMS/lag(PAYEMS,1)-1)
Emp <- Emp  %>% mutate (norm_rate = scale(rate, center = TRUE, scale = TRUE))
Emp <- select(Emp, -c(rate))

GDP <- GDP %>% mutate(rate =GDPC1/lag(GDPC1,1)-1)
GDP <- GDP  %>% mutate (norm_rate = scale(rate, center = TRUE, scale = TRUE))
GDP <- select(GDP, -c(rate,GDPC1))

months <- lapply(X = GDP$date, FUN = seq.Date, by = "month", length.out = 3)
months <- data.frame(date = do.call(what = c, months))

m_GDP <- left_join(x = months, y = GDP , by = "date")

df <- cbind(m_GDP,Emp$norm_rate)

names(df) <- c("date","S01_GDP","S02_Emp")

df_marss  <- df %>% gather(key = "serie", value = "value", -date)

df_marss <- df_marss  %>% spread(key=date,value=value)

df_marss$serie <- NULL

df_marss <- as.matrix(df_marss)

# 3) Define State Space matrices

# Matrix Z

Z <- matrix(list("0.33*z1","z2",
                 "0.67*z1",0,
                 "z1",0,
                 "0.67*z1",0,
                 "0.33*z1",0,
                 1/3,0,
                 2/3,0,
                 1,0,
                 2/3,0,
                 1/3,0,
                 0,1,
                 0,0),2,12)

m <- nrow(Z)
p <- ncol(Z)

# Matrix R

R <- matrix(list(0),m,m)

# Matrix B

B <- matrix (list("b1",1,0,0,0,0,0,0,0,0,0,0,
                  "b2",0,1,0,0,0,0,0,0,0,0,0,
                  0,0,0,1,0,0,0,0,0,0,0,0,
                  0,0,0,0,1,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"b6",1,0,0,0,0,0,
                  0,0,0,0,0,"b7",0,1,0,0,0,0,
                  0,0,0,0,0,0,0,0,1,0,0,0,
                  0,0,0,0,0,0,0,0,0,1,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"b11",1,
                  0,0,0,0,0,0,0,0,0,0,"b12",0),12,12) 


# Matrix Q

Q <- matrix (list("q1",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"q6",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"q11",0,
                  0,0,0,0,0,0,0,0,0,0,0,0),12,12)

# Other options for matrix Q

# Q <- matrix(list(0),12,12)
# Q <- ldiag(list("q1", 0,0,0,0,"q6",0,0,0,0,"q11",0))

# Q <- "diagonal and unequal"

# Rest of matrices

x0 <-  matrix(0,p,1)
A  <- matrix(0,(length(df)-1),1)
U <- matrix(0,p,1)
V0 <- 5*diag(1,p)
U <-  matrix(0,p,1)


# 4) Estimation

# Define model

model.gen =list(Z=Z,A=A,R=R,B=B,U=U,Q=Q,x0=x0,V0=V0,tinitx=0)

# Estimation

kf_ss= try(MARSS(df_marss, model=model.gen,control= list(trace=1,maxit = 300), method="BFGS"), silent=TRUE)

test_that("GDF example for numerical stabilty", {
  expect_true(!inherits(kf_ss, "try-error"))
})

kf_ss= try(MARSS(df_marss, model=model.gen,control= list(trace=0, maxit = 300), method="BFGS"), silent=TRUE)

