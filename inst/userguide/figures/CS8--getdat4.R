###################################################
### code chunk number 23: getdat4
###################################################
years <- rockfish[, 1]
dat <- rockfish[, -1]
dat4 <- dat[years >= 1980, ]
years4 <- years[years >= 1980]
isdata <- apply(is.na(dat4), 2, sum) != length(years4)
dat4 <- dat4[, isdata]


