###################################################
### code chunk number 15: paramvals
###################################################
den91.params <- matrix(NA, nrow = 11, ncol = 2, dimnames = list(c(paste("sim", 1:9), "mean sim", "true"), c("den91.U", "den91.Q")))
for (i in 1:9) {
  den.years <- years[!is.na(y.tss[i, ])] # the non missing years
  den.yts <- y.tss[i, !is.na(y.tss[i, ])] # the non missing counts
  den.n.yts <- length(den.years)
  delta.pop <- rep(NA, den.n.yts - 1) # create a vector to store transitions
  tau <- rep(NA, den.n.yts - 1) # create a vector of time step sizes
  for (t in 2:den.n.yts) {
    delta.pop[t - 1] <- den.yts[t] - den.yts[t - 1] # store each transition
    tau[t - 1] <- den.years[t] - den.years[t - 1] # the transitions
  } # end t loop
  den91 <- lm(delta.pop ~ -1 + tau) # note: the "-1" specifies no intercept
  den91.params[i, ] <- c(den91$coefficients, var(resid(den91)))
}
den91.params[10, ] <- apply(den91.params[1:9, ], 2, mean)
den91.params[11, ] <- c(sim.u, sim.Q)


