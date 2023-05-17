###################################################
### code chunk number 12: Covar_sec6_02_set-up-month-factors
###################################################
# number of "seasons" (e.g., 12 months per year)
period <- 12
# first "season" (e.g., Jan = 1, July = 7)
per.1st <- 1
# create factors for seasons
c.in <- diag(period)
for (i in 2:(ceiling(TT / period))) {
  c.in <- cbind(c.in, diag(period))
}
# trim c.in to correct start & length
c.in <- c.in[, (1:TT) + (per.1st - 1)]
# better row names
rownames(c.in) <- month.abb


