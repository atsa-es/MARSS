###################################################
### code chunk number 4: Cs2-readindata
###################################################
dat <- t(harborSealWA) # Transpose
years <- dat[1, ] # [1,] means row 1
n <- nrow(dat) - 1
dat <- dat[2:nrow(dat), ] # no years


