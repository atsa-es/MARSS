###################################################
### code chunk number 7: Cs1_Exercise1
###################################################
par(mfrow = c(3, 3))
sim.u <- -0.05
sim.Q <- 0.02
sim.R <- 0.05
nYr <- 50
fracmiss <- 0.1
init <- 7
years <- seq(1:nYr)
for (i in 1:9) {
  x <- rep(NA, nYr) # vector for ts w/o measurement error
  y <- rep(NA, nYr) # vector for ts w/ measurement error
  x[1] <- init
  for (t in 2:nYr) {
    x[t] <- x[t - 1] + sim.u + rnorm(1, mean = 0, sd = sqrt(sim.Q))
  }
  for (t in 1:nYr) {
    y[t] <- x[t] + rnorm(1, mean = 0, sd = sqrt(sim.R))
  }
  missYears <-
    sample(years[2:(nYr - 1)], floor(fracmiss * nYr), replace = FALSE)
  y[missYears] <- NA
  plot(years, y,
    xlab = "", ylab = "log abundance", lwd = 2, bty = "l"
  )
  lines(years, x, type = "l", lwd = 2, lty = 2)
  title(paste("simulation ", i))
}
legend("topright", c("Observed", "True"),
  lty = c(-1, 2), pch = c(1, -1)
)


###################################################
### code chunk number 17: Cs1_Exercise2
###################################################
sim.u <- -0.05 # growth rate
sim.Q <- 0.02 # process error variance
sim.R <- 0.05 # non-process error variance
nYr <- 50 # number of years of data to generate
fracmiss <- 0.1 # fraction of years that are missing
init <- 7 # log of initial pop abundance (~1100 individuals)
nsim <- 9
years <- seq(1:nYr) # col of years
params <- matrix(NA,
  nrow = (nsim + 2), ncol = 5,
  dimnames = list(
    c(paste("sim", 1:nsim), "mean sim", "true"),
    c("kem.U", "den91.U", "kem.Q", "kem.R", "den91.Q")
  )
)
x.ts <- matrix(NA, nrow = nsim, ncol = nYr) # ts w/o measurement error
y.ts <- matrix(NA, nrow = nsim, ncol = nYr) # ts w/ measurement error
for (i in 1:nsim) {
  x.ts[i, 1] <- init
  for (t in 2:nYr) {
    x.ts[i, t] <- x.ts[i, t - 1] + sim.u + rnorm(1, mean = 0, sd = sqrt(sim.Q))
  }
  for (t in 1:nYr) {
    y.ts[i, t] <- x.ts[i, t] + rnorm(1, mean = 0, sd = sqrt(sim.R))
  }
  missYears <- sample(years[2:(nYr - 1)], floor(fracmiss * nYr),
    replace = FALSE
  )
  y.ts[i, missYears] <- NA

  # MARSS estimates
  kem <- MARSS(y.ts[i, ], silent = TRUE)
  # type=vector outputs the estimates as a vector instead of a list
  params[i, c(1, 3, 4)] <- coef(kem, type = "vector")[c(2, 3, 1)]

  # Dennis et al 1991 estimates
  den.years <- years[!is.na(y.ts[i, ])] # the non missing years
  den.yts <- y.ts[i, !is.na(y.ts[i, ])] # the non missing counts
  den.n.yts <- length(den.years)
  delta.pop <- rep(NA, den.n.yts - 1) # transitions
  tau <- rep(NA, den.n.yts - 1) # time step lengths
  for (t in 2:den.n.yts) {
    delta.pop[t - 1] <- den.yts[t] - den.yts[t - 1] # transitions
    tau[t - 1] <- den.years[t] - den.years[t - 1] # time step length
  } # end i loop
  den91 <- lm(delta.pop ~ -1 + tau) # -1 specifies no intercept
  params[i, c(2, 5)] <- c(den91$coefficients, var(resid(den91)))
}
params[nsim + 1, ] <- apply(params[1:nsim, ], 2, mean)
params[nsim + 2, ] <- c(sim.u, sim.u, sim.Q, sim.R, sim.Q)


###################################################
### code chunk number 21: Cs1_Exercise3
###################################################
# Needs Example 2 to be run first
par(mfrow = c(3, 3))
pd <- 0.1
xd <- -log(pd) # decline threshold
te <- 100
tyrs <- 1:te # extinction time horizon
for (j in c(10, 1:8)) {
  real.ex <- denn.ex <- kal.ex <- matrix(nrow = te)

  # MARSS parameter estimates
  u <- params[j, 1]
  Q <- params[j, 3]
  if (Q == 0) Q <- 1e-4 # just so the extinction calc doesn't choke
  p.ever <- ifelse(u <= 0, 1, exp(-2 * u * xd / Q))
  for (i in 1:100) {
    if (is.finite(exp(2 * xd * abs(u) / Q))) {
      sec.part <- exp(2 * xd * abs(u) / Q) * pnorm((-xd - abs(u) * tyrs[i]) / sqrt(Q * tyrs[i]))
    } else {
      sec.part <- 0
    }
    kal.ex[i] <- p.ever * pnorm((-xd + abs(u) * tyrs[i]) / sqrt(Q * tyrs[i])) + sec.part
  } # end i loop

  # Dennis et al 1991 parameter estimates
  u <- params[j, 2]
  Q <- params[j, 5]
  p.ever <- ifelse(u <= 0, 1, exp(-2 * u * xd / Q))
  for (i in 1:100) {
    denn.ex[i] <- p.ever * pnorm((-xd + abs(u) * tyrs[i]) / sqrt(Q * tyrs[i])) +
      exp(2 * xd * abs(u) / Q) * pnorm((-xd - abs(u) * tyrs[i]) / sqrt(Q * tyrs[i]))
  } # end i loop

  # True parameter values
  u <- sim.u
  Q <- sim.Q
  p.ever <- ifelse(u <= 0, 1, exp(-2 * u * xd / Q))
  for (i in 1:100) {
    real.ex[i] <- p.ever * pnorm((-xd + abs(u) * tyrs[i]) / sqrt(Q * tyrs[i])) +
      exp(2 * xd * abs(u) / Q) * pnorm((-xd - abs(u) * tyrs[i]) / sqrt(Q * tyrs[i]))
  } # end i loop

  # plot it
  plot(tyrs, real.ex,
    xlab = "time steps into future",
    ylab = "probability of extinction", ylim = c(0, 1), bty = "l"
  )
  if (j <= 8) title(paste("simulation ", j))
  if (j == 10) title("average over sims")
  lines(tyrs, denn.ex, type = "l", col = "red", lwd = 2, lty = 1)
  lines(tyrs, kal.ex, type = "l", col = "green", lwd = 2, lty = 2)
}
legend("bottomright", c("True", "Dennis", "KalmanEM"),
  pch = c(1, -1, -1),
  col = c(1, 2, 3), lty = c(-1, 1, 2), lwd = c(-1, 2, 2), bty = "n"
)


###################################################
### code chunk number 23: Cs1_Exercise4
###################################################
par(mfrow = c(1, 1))
CSEGtmufigure(N = 50, u = -0.05, s2p = 0.02)


###################################################
### code chunk number 27: Cs1_Exercise5
###################################################
# If you have your data in a tab delimited file with a header
# This is how you would read it in using file.choose()
# to call up a directory browser.
# However, the package has the datasets for the examples
# dat=read.table(file.choose(), skip=1)
# dat=as.matrix(dat)
dat <- wilddogs
CSEGriskfigure(dat, CI.method = "hessian", silent = TRUE)


###################################################
### code chunk number 2: Cs1_a_required_libraries
###################################################
library(MARSS)


