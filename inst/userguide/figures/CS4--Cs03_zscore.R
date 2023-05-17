###################################################
### code chunk number 5: Cs03_zscore
###################################################
Sigma <- sqrt(apply(dat.spp.1980, 1, var, na.rm = TRUE))
y.bar <- apply(dat.spp.1980, 1, mean, na.rm = TRUE)
dat.z <- (dat.spp.1980 - y.bar) * (1 / Sigma)
rownames(dat.z) <- rownames(dat.spp.1980)


