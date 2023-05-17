###################################################
### code chunk number 4: Cs02_transpose_data
###################################################
# transpose data so time goes across columns
dat.spp.1980 <- t(dat.spp.1980)
N.ts <- nrow(dat.spp.1980)
TT <- ncol(dat.spp.1980)


