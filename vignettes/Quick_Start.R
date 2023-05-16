## -----------------------------------------------------------------------------
B1 <- matrix(list("b",0,0,"b"),2,2)
U1 <- matrix(0,2,1)
Q1 <- matrix(c("q11","q12","q12","q22"),2,2)
Z1 <- matrix(c(1,0,1,1,1,0),3,2)
A1 <- matrix(list("a1",0,0),3,1)
R1 <- matrix(list("r11",0,0,0,"r",0,0,0,"r"),3,3)
pi1 <- matrix(0,2,1); V1=diag(1,2)
model.list <- list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)

## ----eval=FALSE---------------------------------------------------------------
#  fit <- MARSS(y, model=model.list)

## -----------------------------------------------------------------------------
library(MARSS)
set.seed(1234)
x <- rbind(arima.sim(n=50,list(ar=0.95), sd=0.4), 
           arima.sim(n=50,list(ar=0.95), sd=.02))
y <- Z1 %*% x + matrix(rnorm(3*50,0,0.1), 3, 50)
fit <- MARSS(y, model=model.list, silent=TRUE)
tidy(fit)

## -----------------------------------------------------------------------------
fit1 <- MARSS(y, model=model.list, silent=TRUE)
fit2 <- MARSS(y, model=model.list, method="BFGS", silent=TRUE)
fit3 <- MARSS(y, model=model.list, method="TMB", silent=TRUE)
cat("logLL of the 3 fits:", c(fit1$logLik, fit2$logLik, fit3$logLik))

## ----results="hide"-----------------------------------------------------------
fit1 <- MARSS(y, model=model.list, control = list(maxit=15))
fit2 <- MARSS(y, model=model.list, method="BFGS", inits = fit1)

## -----------------------------------------------------------------------------
temp <- matrix(rnorm(50, seq(0,1,1/50), 0.1),nrow=1)
C1 <- matrix(c("temp1","temp2"),2,1)
model.list$C <- C1
model.list$c <- temp

## ----results="hide"-----------------------------------------------------------
fit <- MARSS(y, model=model.list, method="BFGS")

## -----------------------------------------------------------------------------
yts <- ts(t(y), frequency = 12) # requires time down the rows
fcov <- forecast::fourier(yts,1) |> t()

## -----------------------------------------------------------------------------
yts <- ts(t(y), frequency = 12) # monthly data
mcov <- forecast::seasonaldummy(yts) |> t() # month factor

## -----------------------------------------------------------------------------
B1 <- matrix(list("-0.1+1*b",0,0,"0.1+1*b"),2,2)
Q1 <- matrix(list(1,0,0,1),2,2)
Z1 <- matrix(list("1*z1+-1*z2",0,"z2","2*z1","z1",0),3,2)
model.list <- list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=0)

## -----------------------------------------------------------------------------
fit <- MARSS(y, model=model.list, silent = TRUE)

## -----------------------------------------------------------------------------
TT <- dim(y)[2]
B1 <- array(list(),dim=c(2,2,TT))
B1[,,1:20] <- matrix(list("b",0,0,"b_1"),2,2)
B1[,,21:TT] <- matrix(list("b",0,0,"b_2"),2,2)

