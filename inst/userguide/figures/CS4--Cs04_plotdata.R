###################################################
### code chunk number 7: Cs04_plotdata
###################################################
spp <- rownames(dat.spp.1980)
par(mfcol = c(3, 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in spp) {
  plot(dat.z[i, ], xlab = "", ylab = "Abundance index", bty = "L", xaxt = "n", pch = 16, col = "blue", type = "b")
  axis(1, 12 * (0:dim(dat.spp.1980)[2]) + 1, 1980 + 0:dim(dat.spp.1980)[2])
  title(i)
}


