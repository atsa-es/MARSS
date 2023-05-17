###################################################
### code chunk number 9: Cs01_set.up.data
###################################################
years <- harborSeal[, 1] # first col is years
# leave off Hood Canal data for now
sealData <- t(harborSeal[, c(2:7, 9:13)])


