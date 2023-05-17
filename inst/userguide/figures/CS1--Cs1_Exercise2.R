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


