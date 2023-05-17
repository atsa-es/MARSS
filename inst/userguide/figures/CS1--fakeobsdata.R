###################################################
### code chunk number 5: fakeobsdata
###################################################
for (t in 1:nYr) {
  y[t] <- x[t] + rnorm(1, mean = 0, sd = sqrt(sim.R))
}
missYears <- sample(years[2:(nYr - 1)], floor(fracmissing * nYr),
  replace = FALSE
)
y[missYears] <- NA


