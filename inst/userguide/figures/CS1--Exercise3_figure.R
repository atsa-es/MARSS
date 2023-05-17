###################################################
### code chunk number 20: Exercise3_figure
###################################################
# Needs Example 2 to be run first
par(mfrow = c(3, 3))
pd <- 0.1; xd <- -log(pd)
te <- 100; tyrs <- 1:te # extinction time horizon
for (j in c(10, 1:8)) {
  real.ex <- matrix(nrow = te)
  denn.ex <- matrix(nrow = te)
  kal.ex <- matrix(nrow = te)

  # MARSS
  u <- params[j, 1]; Q <- params[j, 3]
  if (Q == 0) Q <- 1e-4 # just so the extinction calc doesn't choke
  p.ever <- ifelse(u <= 0, 1, exp(-2 * u * xd / Q))
  for (i in 1:100) {
    if (is.finite(exp(2 * xd * abs(u) / Q))) { # Q!=0
      part2 <- exp(2 * xd * abs(u) / Q) * 
               pnorm((-xd - abs(u) * tyrs[i]) / sqrt(Q * tyrs[i]))
    } else { part2 <- 0 } # Q=0
    kal.ex[i] <- p.ever * 
                 pnorm((-xd + abs(u) * tyrs[i]) / sqrt(Q * tyrs[i])) + 
                 part2
  } # end i loop

  # Dennis et al 1991
  u <- params[j, 2]; Q <- params[j, 5]
  p.ever <- ifelse(u <= 0, 1, exp(-2 * u * xd / Q))
  for (i in 1:100) {
    denn.ex[i] <- p.ever * pnorm((-xd + abs(u) * tyrs[i]) / (sqrt(Q) * sqrt(tyrs[i]))) + 
      exp(2 * xd * abs(u) / Q) * pnorm((-xd - abs(u) * tyrs[i]) / (sqrt(Q) * sqrt(tyrs[i])))
  } # end i loop

  # True
  u <- sim.u; Q <- sim.Q
  p.ever <- ifelse(u <= 0, 1, exp(-2 * u * xd / Q))
  for (i in 1:100) {
    real.ex[i] <- p.ever * pnorm((-xd + abs(u) * tyrs[i]) / sqrt(Q * tyrs[i])) + 
      exp(2 * xd * abs(u) / Q) * pnorm((-xd - abs(u) * tyrs[i]) / sqrt(Q * tyrs[i]))
  } # end i loop

  plot(tyrs, real.ex, xlab = "Time steps into future", 
       ylab = "Probability of extinction", ylim = c(0, 1), bty = "l")
  if (j <= 8) title(paste("simulation ", j))
  if (j == 10) title("average over sims")
  lines(tyrs, denn.ex, type = "l", col = "red", lwd = 2, lty = 1) 
  lines(tyrs, kal.ex, type = "l", col = "green", lwd = 2, lty = 2)
}
legend("bottomright", c("True", "Dennis", "KalmanEM"), pch = c(1, -1, -1), 
       col = c(1, 2, 3), lty = c(-1, 1, 2), lwd = c(-1, 2, 2), bty = "n")


