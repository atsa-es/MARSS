###################################################
### code chunk number 2: Cs00_required-libraries
###################################################
library(MARSS)
library(tidyr)
library(ggplot2)
library(forecast)


###################################################
### code chunk number 3: Cs101_structTS-level
###################################################
y <- window(treering, start = 0, end = 20)

fit1 <- StructTS(y, type = "level")


###################################################
### code chunk number 4: Cs102_structTS-level
###################################################
vy <- var(y, na.rm = TRUE) / 100
mod.list <- list(
  x0 = matrix(y[1]), U = "zero", tinitx = 0,
  Q = matrix(fit1$coef[1]), R = matrix(fit1$coef[2]),
  V0 = matrix(1e+06 * vy)
)
fit2 <- MARSS(as.vector(y), model = mod.list)
# Now estimate the parameters
mod.list <- list(
  x0 = matrix(y[1]), U = "zero", tinitx = 0, V0 = matrix(1e+06 * vy),
  Q = matrix("s2xi"), R = matrix("s2eps")
)
fit3 <- MARSS(as.vector(y), model = mod.list, method = "BFGS")
fit4 <- MARSS(as.vector(y), model = mod.list, 
              control = list(allow.degen = FALSE))


###################################################
### code chunk number 5: Cs103_structTS-level
###################################################
fit2$kf <- MARSSkfss(fit2)
fit3$kf <- MARSSkfss(fit3)
fit4$kf <- MARSSkfss(fit4)
df <- data.frame(
  StructTS = fit1$fitted, fit2 = fit2$kf$xtt[1, ],
  fit.bfgs = fit3$kf$xtt[1, ], fit.em = fit4$kf$xtt[1, ]
)
head(df)


###################################################
### code chunk number 6: Cs104_structTS-level
###################################################
require(tidyr)
require(ggplot2)
df1 <- as.data.frame(fit1$fitted)
vars <- colnames(df1)
df2 <- as.data.frame(t(fit3$kf$xtt))
colnames(df2) <- vars
df3 <- as.data.frame(t(fit4$kf$xtt))
colnames(df3) <- vars
df1$model <- "StructTS"
df2$model <- "MARSS BFGS"
df3$model <- "MARSS EM"
df1$t <- as.vector(time(fit1$fitted))
df2$t <- df1$t
df3$t <- df1$t
df <- rbind(df1, df2, df3)
df <- df %>% pivot_longer(all_of(vars))
ggplot(df, aes(x = t, y = value)) +
  geom_line() +
  facet_wrap(~model) +
  ggtitle("Level estimate from model fit with StructTS and MARSS")


###################################################
### code chunk number 7: Cs201_structTS-leveltrend
###################################################
y <- log10(forecast:::subset.ts(UKgas, quarter = 2))
fit1 <- StructTS(y, type = "trend")


###################################################
### code chunk number 8: Cs202_structTS-leveltrend
###################################################
vy <- var(y, na.rm = TRUE) / 100
B <- matrix(c(1, 0, 1, 1), 2, 2)
Z <- matrix(c(1, 0), 1, 2)
# fitx parameters at fit1 values
mod.list <- list(
  x0 = matrix(c(y[1], 0), 2, 1), U = "zero", tinitx = 0,
  Q = diag(fit1$coef[1:2]), R = matrix(fit1$coef[3]),
  V0 = matrix(1e+06 * vy, 2, 2), Z = Z, B = B
)
fit2 <- MARSS(as.vector(y), model = mod.list, fit = FALSE, 
              control = list(trace = -1))
fit2$par <- fit2$start # otherwise par is NULL since fit=FALSE


###################################################
### code chunk number 9: Cs203_structTS-leveltrend
###################################################
mod.list <- list(
  x0 = matrix(c(y[1], 0), 2, 1), U = "zero", tinitx = 0,
  Q = ldiag(c("s2xi", "s2zeta")), R = matrix("s2eps"),
  V0 = matrix(1e+06 * vy, 2, 2) + diag(1e-10, 2), Z = Z, B = B
)
fit3 <- MARSS(as.vector(y), model = mod.list, method = "BFGS")
fit4 <- MARSS(as.vector(y), model = mod.list, 
              control = list(allow.degen = FALSE))


###################################################
### code chunk number 10: Cs204_structTS-leveltrend
###################################################
fit2$kf <- MARSSkfss(fit2)
fit3$kf <- MARSSkfss(fit3)
fit4$kf <- MARSSkfss(fit4)
data.frame(
  StructTS = fit1$fitted[, 2], fit2 = fit2$kf$xtt[2, ],
  fit.bfgs = fit3$kf$xtt[2, ], fit.em = fit4$kf$xtt[2, ]
)[1:5, ]


###################################################
### code chunk number 11: Cs205_structTS-leveltrend
###################################################
require(tidyr)
require(ggplot2)
df1 <- as.data.frame(fit1$fitted)
vars <- colnames(df1)
df2 <- as.data.frame(t(fit3$kf$xtt))
colnames(df2) <- vars
df3 <- as.data.frame(t(fit4$kf$xtt))
colnames(df3) <- vars
df1$model <- "StructTS"
df2$model <- "MARSS BFGS"
df3$model <- "MARSS EM"
df1$t <- as.vector(time(fit1$fitted))
df2$t <- df1$t
df3$t <- df1$t
df <- rbind(df1, df2, df3)
df <- df %>% pivot_longer(all_of(vars))
ggplot(df, aes(x = t, y = value, color = model, linetype = model, shape = model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free") +
  scale_linetype_manual("model", values = c(1, 1, 0)) +
  scale_shape_manual("model", values = c(NA, NA, 16))


###################################################
### code chunk number 12: Cs301_structTS-bsm
###################################################
y <- log10(UKgas)
fit1 <- StructTS(y, type = "BSM")


###################################################
### code chunk number 13: Cs302_structTS-bsm
###################################################
makeB <- function(nf) {
  B <- matrix(0, nf + 1L, nf + 1L)
  B[1L:2L, 1L:2L] <- c(1, 0, 1, 1)
  B[3L, ] <- c(0, 0, rep(-1, nf - 1L))
  if (nf >= 3L) {
    ind <- 3:nf
    B[cbind(ind + 1L, ind)] <- 1
  }
  return(B)
}


###################################################
### code chunk number 14: Cs303_structTS-bsm
###################################################
nf <- frequency(y)
vy <- var(y) / 100
B <- makeB(nf)
Z <- matrix(c(1, 0, 1, rep(0, nf - 2L)), 1, nf + 1)

Q <- ldiag(list("s2xi", "s2zeta", "s2w", 0, 0))
R <- matrix("s2eps")
V0 <- matrix(1e+06 * vy, nf + 1, nf + 1) + diag(1e-10, nf + 1)
mod.list <- list(
  x0 = matrix(c(y[1], rep(0, nf)), ncol = 1), 
  U = "zero", A = "zero", tinitx = 0,
  Q = Q, R = R, V0 = V0, Z = Z, B = B
)
fit3 <- MARSS(as.vector(y), model = mod.list, method = "BFGS")
fit4 <- MARSS(as.vector(y), model = mod.list, 
              control = list(allow.degen = FALSE))
fit4$kf <- MARSSkfss(fit4)
fit3$kf <- MARSSkfss(fit3)


###################################################
### code chunk number 15: Cs304_structTS-bsm
###################################################
require(tidyr)
require(ggplot2)
df1 <- as.data.frame(fit1$fitted)
vars <- colnames(df1)
df2 <- as.data.frame(t(fit3$kf$xtt)[, 1:3])
colnames(df2) <- vars
df3 <- as.data.frame(t(fit4$kf$xtt)[, 1:3])
colnames(df3) <- vars
df1$model <- "StructTS"
df2$model <- "MARSS BFGS"
df3$model <- "MARSS EM"
df1$t <- as.vector(time(fit1$fitted))
df1$Qtr <- as.vector(cycle(fit1$fitted))
df2$t <- df1$t
df2$Qtr <- df1$Qtr
df3$t <- df1$t
df3$Qtr <- df1$Qtr
df <- rbind(df1, df2, df3)
df <- subset(df, Qtr == 1) %>% pivot_longer(all_of(vars))
ggplot(df, aes(x = t, y = value, color = model, linetype = model, shape = model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free") +
  scale_linetype_manual("model", values = c(1, 1, 0)) +
  scale_shape_manual("model", values = c(NA, NA, 15))


###################################################
### code chunk number 16: Cs401_forecast
###################################################
y <- log10(UKgas)
fit1 <- StructTS(y, type = "BSM")

nf <- frequency(y)
vy <- var(y) / 100
B <- makeB(nf) # defined in the BSM section above
Z <- matrix(c(1, 0, 1, rep(0, nf - 2L)), 1, nf + 1)
V0 <- matrix(1e+06 * vy, nf + 1, nf + 1) + diag(1e-10, nf + 1)
mod.list <- list(
  x0 = matrix(c(y[1], rep(0, nf)), ncol = 1), U = "zero", A = "zero", tinitx = 0,
  Q = diag(c(fit1$coef[1:3], 0, 0)), R = matrix(fit1$coef[4]), V0 = V0, Z = Z, B = B
)
fit2 <- MARSS(as.vector(y), model = mod.list)


###################################################
### code chunk number 17: Cs402_forecast
###################################################
fr1 <- predict(fit1, n.ahead = 5)
fr1


###################################################
### code chunk number 18: Cs403_forecast
###################################################
fr2 <- predict(fit2, n.ahead = 5, interval="prediction")
fr2


###################################################
### code chunk number 19: Cs404_forecast
###################################################
rbind(
  pred1 = fr1$pred, pred2 = fr2$pred$estimate[fr2$ft],
  se1 = fr1$se, se2 = fr2$pred$se[fr2$ft]
)


###################################################
### code chunk number 20: Cs405_forecast
###################################################
fr1 <- forecast::forecast(fit1, h = 10)
fr2 <- forecast::forecast(fit2, h = 10)
p1 <- ggplot2::autoplot(fr1, include = 8)
p2 <- ggplot2::autoplot(fr2, include = 8)
gridExtra::grid.arrange(p1, p2, nrow = 1)


###################################################
### code chunk number 21: Cs501_fitted
###################################################
fitted1 <- fitted(fit1)
plot(fitted1)


###################################################
### code chunk number 22: Cs502_fitted
###################################################
fitted2 <- tsSmooth(fit2, type = "xtt")
fitted2 <- subset(fitted2, .rownames  %in% c("X1", "X2", "X3"))


###################################################
### code chunk number 23: Cs503_fitted
###################################################
ggplot(fitted2, aes(x = t, y = .estimate)) +
  geom_line() +
  facet_wrap(~.rownames, ncol=1, scale="free_y")


###################################################
### code chunk number 24: Cs504_fitted
###################################################
fitted3 <- MARSSkfss(fit2)$xtt
fitted3 <- ts(t(fitted3[1:3,]))
plot(fitted3)


###################################################
### code chunk number 25: Cs505_fitted
###################################################
fitted2 <- fitted(fit2, type="ytt")
df2 <- data.frame(t=as.numeric(time(fitted1)), fitted=fitted2$.fitted, name="MARSS")
df1 <- data.frame(t=df2$t, fitted=as.numeric(fitted1[,1]+fitted1[,3]), name="StructTS")
df <- rbind(df1, df2)
df$y <- fitted2$y

ggplot(df) +
  geom_line(aes(x = t, y = fitted)) + geom_point(aes(x = t, y = y), col="blue") +
  facet_wrap(~name, ncol=1)


###################################################
### code chunk number 26: Cs601_residuals
###################################################
resids1 <- residuals(fit1)


###################################################
### code chunk number 27: Cs602_residuals
###################################################
resids2 <- residuals(fit2, type = "tt", standardization="marginal")


###################################################
### code chunk number 28: Cs603_residuals
###################################################
df2 <- data.frame(t=as.numeric(time(resids1)), resids=resids2$.std.resids, name="MARSS")
df1 <- data.frame(t=df2$t, resids=as.numeric(resids1), name="StructTS")
df3 <- data.frame(t=df2$t, resids=df1$resids-df2$resids, name="difference")
df <- rbind(df1, df2, df3)

ggplot(df, aes(x = t, y = resids)) +
  geom_line() + facet_wrap(~name, ncol=1, scale="free_y") +
  ggtitle("Marginal standardized model residuals")


###################################################
### code chunk number 29: Cs701_multivariate
###################################################
set.seed(100)
TT <- 60; t <- 1:TT
q <- 0.01; r <- 0.01
trend <- 0.2*sin((1:TT)/4)
level <- cumsum(rnorm(TT,trend, sqrt(q)))

# Simulated data
n <- 5
miss.percent <- 0.5
ym <- matrix(1,n,1)%*%level + matrix(rnorm(TT*n,0,sqrt(r*100)),n,TT)
ym[sample(n*TT, miss.percent*n*TT)] <- NA


###################################################
### code chunk number 30: Cs702_multivariate
###################################################
par(mfrow=c(2,1), mar=c(3,3,1,1))
ylims <- c(min(ym, na.rm=TRUE),max(ym, na.rm=TRUE))
plot(t, trend, ylim=ylims, col="red", type="l"); lines(t, level, col="black")
legend("topright", c("trend", "level"), lty=1, col=c("red", "black"))
matplot(t,t(ym), pch=1:n, ylab="y", xlab="", ylim=ylims, main="bad data")
lines(t,level)


###################################################
### code chunk number 31: Cs703_multivariate
###################################################
vy <- var(y, na.rm = TRUE) / 100
mod.list.x <- list(
  x0 = matrix(list("x0",0),nrow=2), tinitx = 1,
  V0 = matrix(1e+06 * vy, 2, 2) + diag(1e-10, 2),
  Q = ldiag(list(q,"qt")),
  B = matrix(c(1, 0, 1, 1), 2, 2),
  U = "zero")


###################################################
### code chunk number 32: Cs704_multivariate
###################################################
mod.list.y <- list(
  A = "zero",
  R = "diagonal and equal")


###################################################
### code chunk number 33: Cs705_multivariate
###################################################
Z <- matrix(c(1, 0), 1, 2, byrow=TRUE)
mod.list <- c(mod.list.x, mod.list.y, list(Z=Z))
fitu <- MARSS(ym[1,], model = mod.list, method = "BFGS", inits=list(x0=0))


###################################################
### code chunk number 34: Cs706_multivariate
###################################################
Z <- matrix(c(1, 0), n, 2, byrow=TRUE)
mod.list <- c(mod.list.x, mod.list.y, list(Z=Z))
fitm <- MARSS(ym, model = mod.list, method = "BFGS", inits=list(x0=0))


###################################################
### code chunk number 35: Cs707_multivariate
###################################################
true <- data.frame(.rownames=rep(c("X1","X2"), each=TT), t=t, 
                   .estimate=c(level,trend), .se=NA, name="true")
statesu <- tsSmooth(fitu)
statesu$name <- "one bad"
statesm <- tsSmooth(fitm)
statesm$name <- "multiple bad"
df <- rbind(true, statesu, statesm)

ggplot(df, aes(x=t, y=.estimate, col=name)) + 
  geom_line() +
  facet_wrap(~.rownames, scale="free_y")


###################################################
### code chunk number 36: Cs708_multivariate
###################################################
par(mfrow=c(2,1), mar=c(3,3,1,1))
covariate <- matrix(c(rep(0,TT-10),rep(1,10)), nrow=1)
ymc <- ym
D <- matrix(c(-1,-1,0,1,1),ncol=1)
ymc <- ym + D%*%covariate
matplot(t,t(ymc), pch=1:n, ylab="y", xlab="", main="data")
lines(level)
plot(t, covariate[1,], col="blue", lty=2, type="l", main="covariate")


###################################################
### code chunk number 37: Cs709_multivariate
###################################################
Z <- matrix(c(1, 0), n, 2, byrow=TRUE)
mod.list <- c(mod.list.x, mod.list.y, list(Z=Z, d=covariate))
fitmc <- MARSS(ymc, model = mod.list, method = "BFGS", inits=list(x0=0))


###################################################
### code chunk number 38: Cs710_multivariate
###################################################
dvals <- data.frame(x=paste0("y",1:n), 
                    val=c(D, coef(fitmc, type="matrix")$D),
                    name=rep(c("true","estimate"), each=n))
ggplot(dvals, aes(x=x, y=val, col=name)) + geom_point() +
  xlab("observation series") + ylab("D estimate") +
  ggtitle("D true and estimated values")


###################################################
### code chunk number 39: Cs711_multivariate
###################################################
statesmc <- tsSmooth(fitmc)
statesmc$name <- "multiple w covariate"
df <- rbind(true, statesu, statesm, statesmc)

ggplot(df, aes(x=t, y=.estimate, col=name)) + 
  geom_line() +
  facet_wrap(~.rownames, scale="free_y")


###################################################
### code chunk number 40: Cs712_multivariate
###################################################
r2 <- r*c(100,10,10,200,400)
a <- runif(n,-1,1)
err <- rnorm(n*TT, mean=rep(a, each=TT), sd=rep(sqrt(r2), each=TT))
ym2 <- matrix(1,nrow=n)%*%level + matrix(err,nrow=n,byrow=TRUE)
ym2[sample(n*TT, miss.percent*n*TT)] <- NA
matplot(t,t(ym2), pch=1:n, ylab="y", xlab="", main="data with different error and bias")
lines(level)


###################################################
### code chunk number 41: Cs713_multivariate
###################################################
Z <- matrix(c(1, 0), n, 2, byrow=TRUE)
mod.list <- c(mod.list.x, list(Z=Z, R="diagonal and unequal", A="scaling"))
fitm2 <- MARSS(ym2, model = mod.list, method = "BFGS", inits=list(x0=0))


###################################################
### code chunk number 42: Cs714_multivariate
###################################################
rvals <- data.frame(x=paste0("y",1:n), 
                    val=c(r2, coef(fitm2)$R),
                    name=rep(c("true","estimate"), each=n))
ggplot(rvals, aes(x=x, y=val, col=name)) + geom_point() +
  xlab("observation series") + ylab("R variance estimate") +
  ggtitle("R true and estimated values")


###################################################
### code chunk number 43: Cs715_multivariate
###################################################
statesm2 <- tsSmooth(fitm2)
statesm2$name <- "multiple w different Rs"
df <- rbind(true, statesu, statesm, statesmc, statesm2)

ggplot(df, aes(x=t, y=.estimate, col=name)) + 
  geom_line() +
  facet_wrap(~.rownames, scale="free_y")


###################################################
### code chunk number 44: Cs716_multivariate
###################################################
set.seed(100)
TT <- 60; t <- 1:TT
q <- 0.5; qt <- 0.01; r <- 0.1
b <- 0.5
trend <- 0.2*sin((1:TT)/4)
level1 <- cumsum(rnorm(TT,trend, sqrt(q)))
level2 <- cumsum(rnorm(TT,trend, sqrt(q)))

# Simulated data
ym <- rbind(level1, level2) + matrix(rnorm(TT*2,0,sqrt(r)),2,TT)


###################################################
### code chunk number 45: Cs717_multivariate
###################################################
par(mfrow=c(2,1), mar=c(3,3,1,1))
ylims <- c(min(ym, na.rm=TRUE),max(ym, na.rm=TRUE))
plot(t, ym[1,], ylim=ylims, type="p"); lines(t, level1, col="black")
plot(t, ym[2,], ylim=ylims, type="p"); lines(t, level2, col="black")


###################################################
### code chunk number 46: Cs718_multivariate
###################################################
vy <- var(y, na.rm = TRUE) / 100
Z <- matrix(c(1, 0), 1, 2)
mod.list.x <- list(
  x0 = matrix(list("x0",0),nrow=2), tinitx = 1,
  V0 = matrix(1e+06 * vy, 2, 2) + diag(1e-10, 2),
  Q = ldiag(list(q,"qt")),
  B = matrix(c(1, 0, 1, 1), 2, 2),
  U = "zero")
mod.list <- c(mod.list.x, mod.list.y, list(Z=Z))
fitm1 <- MARSS(ym[1,], model = mod.list, method = "BFGS", inits=list(x0=0))
fitm2 <- MARSS(ym[2,], model = mod.list, method = "BFGS", inits=list(x0=0))


###################################################
### code chunk number 47: Cs719_multivariate
###################################################
Z <- matrix(c(1, 0, 0, 0, 1, 0), 2, 3, byrow=TRUE)
m <- 3
mod.list.x <- list(
  x0 = matrix(list("x0.1","x0.2",0),nrow=m), tinitx = 1,
  V0 = matrix(1e+06 * vy, m, m) + diag(1e-10, m),
  Q = ldiag(list("q","q","qt")),
  B = matrix(c(1, 0, 1, 0, 1, 1, 0, 0, 1), m, m, byrow=TRUE),
  U = "zero")
mod.list <- c(mod.list.x, mod.list.y, list(Z=Z))
fitm3 <- MARSS(ym, model = mod.list, method = "BFGS", inits=list(x0=0))


###################################################
### code chunk number 48: Cs720_multivariate
###################################################
true <- data.frame(.rownames=rep(c("level 1","level 2","trend"), each=TT), t=t, 
                   .estimate=c(level1,level2,trend), .se=NA, name="true")
statesm1 <- tsSmooth(fitm1)
statesm1$name <- "ts 1 alone"
statesm1$.rownames[statesm1$.rownames=="X2"] <- "trend"
statesm1$.rownames[statesm1$.rownames=="X1"] <- "level 1"
statesm2 <- tsSmooth(fitm2)
statesm2$name <- "ts 2 alone"
statesm2$.rownames[statesm2$.rownames=="X2"] <- "trend"
statesm2$.rownames[statesm2$.rownames=="X1"] <- "level 2"
statesm3 <- tsSmooth(fitm3)
statesm3$name <- "ts 1 & 2 together"
statesm3$.rownames[statesm3$.rownames=="X3"] <- "trend"
statesm3$.rownames[statesm3$.rownames=="X1"] <- "level 1"
statesm3$.rownames[statesm3$.rownames=="X2"] <- "level 2"
df <- rbind(true, statesm1, statesm2, statesm3)

ggplot(df, aes(x=t, y=.estimate, col=name)) + 
  geom_line() +
  facet_wrap(~.rownames, scale="free_y", ncol=1)


