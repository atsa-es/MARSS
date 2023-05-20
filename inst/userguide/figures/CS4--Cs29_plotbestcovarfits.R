###################################################
### code chunk number 42: Cs29_plotbestcovarfits
###################################################
par.mat <- coef(kemz.temp, type = "matrix")
fit.b <- par.mat$Z %*% kemz.temp$states + matrix(par.mat$A, nrow = N.ts, ncol = TT)
spp <- rownames(dat.z)
par(mfcol = c(3, 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:length(spp)) {
  plot(dat.z[i, ], xlab = "", ylab = "Abundance index", bty = "L", xaxt = "n", ylim = c(-4, 3), pch = 16, col = "blue")
  axis(1, 12 * (0:dim(dat.z)[2]) + 1, 1980 + 0:dim(dat.z)[2])
  lines(fit.b[i, ], lwd = 2)
  title(spp[i])
}


