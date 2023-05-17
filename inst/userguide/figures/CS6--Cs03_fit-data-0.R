###################################################
### code chunk number 5: Cs03_fit-data-0
###################################################
# The data is in a ts format, and we need a matrix
dat <- t(as.matrix(Nile))
rownames(dat) <- "Nile"

kem.0 <- MARSS(dat, model = mod.nile.0, silent = TRUE)
summary(kem.0)


