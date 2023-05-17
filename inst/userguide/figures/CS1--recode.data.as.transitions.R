###################################################
### code chunk number 12: recode.data.as.transitions
###################################################
den.years <- years[!is.na(y)] # the non missing years
den.y <- y[!is.na(y)] # the non missing counts
den.n.y <- length(den.years)
delta.pop <- rep(NA, den.n.y - 1) # population transitions
tau <- rep(NA, den.n.y - 1) # step sizes
for (i in 2:den.n.y) {
  delta.pop[i - 1] <- den.y[i] - den.y[i - 1]
  tau[i - 1] <- den.years[i] - den.years[i - 1]
} # end i loop


