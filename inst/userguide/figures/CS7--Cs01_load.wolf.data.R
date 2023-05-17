###################################################
### code chunk number 2: Cs01_load.wolf.data
###################################################
yr1960to2011 <- isleRoyal[, "Year"] >= 1960 & isleRoyal[, "Year"] <= 2011
royale.dat <- log(t(isleRoyal[yr1960to2011, c("Wolf", "Moose")]))


