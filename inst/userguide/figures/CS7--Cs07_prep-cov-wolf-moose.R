###################################################
### code chunk number 9: Cs07_prep-cov-wolf-moose
###################################################
clim.variables <- c(
  "jan.feb.ave.temp", "jan.feb.ave.precip",
  "july.sept.ave.temp"
)
yr1959to2010 <- isleRoyal[, "Year"] >= 1959 & isleRoyal[, "Year"] <= 2010
clim.dat <- t(isleRoyal[yr1959to2010, clim.variables])
z.score.clim.dat <- zscore(clim.dat)


