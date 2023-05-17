###################################################
### code chunk number 15: Covar_sec6_05_month-factor-marss-params
###################################################
# Each taxon has unique density-dependence
B <- "diagonal and unequal"
# Independent process errors
Q <- "diagonal and unequal"
# We have demeaned the data & are fitting a mean-reverting model
# by estimating a diagonal B, thus
U <- "zero"
# Each obs time series is associated with only one process
Z <- "identity"
# The data are demeaned & fluctuate around a mean
A <- "zero"
# Observation errors are independent, but they
# have similar variance due to similar collection methods
R <- "diagonal and equal"
# No covariate effects in the obs equation
D <- "zero"
d <- "zero"


