###################################################
### code chunk number 19: probextcalc (eval = FALSE)
###################################################
## pd <- 0.1 # means a 90 percent decline
## tyrs <- 1:100
## xd <- -log(pd)
## p.ever <- ifelse(u <= 0, 1, exp(-2 * u * xd / Q)) # Q=sigma2
## for (i in 1:100) {
##   Pi[i] <- p.ever * pnorm((-xd + abs(u)*tyrs[i])/sqrt(Q*tyrs[i])) +
##     exp(2*xd*abs(u)/Q) * pnorm((-xd - abs(u)*tyrs[i])/sqrt(Q*tyrs[i]))
## }


