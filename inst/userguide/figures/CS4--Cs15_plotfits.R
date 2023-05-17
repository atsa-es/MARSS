###################################################
### code chunk number 20: Cs15_plotfits
###################################################
fit <- kemz.3
spp <- rownames(dat.z)
par(mfcol = c(3, 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:length(spp)) {
  plot(dat.z[i, ], xlab = "", ylab = "Abundance index", bty = "L", xaxt = "n", ylim = c(-4, 3), pch = 16, col = "blue")
  axis(1, 12 * (0:dim(dat.z)[2]) + 1, 1980 + 0:dim(dat.z)[2])
  par.mat <- coef(fit, type = "matrix")
  lines(as.vector(par.mat$Z[i, , drop = FALSE] %*% fit$states + par.mat$A[i, ]), lwd = 2)
  title(spp[i])
}


