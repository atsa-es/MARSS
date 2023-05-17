###################################################
### code chunk number 17: Covar_sec6_07_poly-month-factor
###################################################
# number of "seasons" (e.g., 12 months per year)
period <- 12
# first "season" (e.g., Jan = 1, July = 7)
per.1st <- 1
# order of polynomial
poly.order <- 3
# create polynomials of months
month.cov <- matrix(1, 1, period)
for (i in 1:poly.order) {
  month.cov <- rbind(month.cov, (1:12)^i)
}
# our c matrix is month.cov replicated once for each year
c.m.poly <- matrix(month.cov, poly.order + 1, TT + period, byrow = FALSE)
# trim c.in to correct start & length
c.m.poly <- c.m.poly[, (1:TT) + (per.1st - 1)]

# Everything else remains the same as in the previous example
model.list <- list(
  B = B, U = U, Q = Q, Z = Z, A = A, R = R,
  C = C, c = c.m.poly, D = D, d = d
)
seas.mod.2 <- MARSS(dat, model = model.list, control = list(maxit = 1500))


