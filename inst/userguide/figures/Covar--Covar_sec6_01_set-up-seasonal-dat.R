###################################################
### code chunk number 11: Covar_sec6_01_set-up-seasonal-dat
###################################################
years <- fulldat[, "Year"] >= 1965 & fulldat[, "Year"] < 1975
phytos <- c(
  "Diatoms", "Greens", "Bluegreens",
  "Unicells", "Other.algae"
)
dat <- t(fulldat[years, phytos])

# z.score data again because we changed the mean when we subsampled
dat <- zscore(dat)
# number of time periods/samples
TT <- ncol(dat)


