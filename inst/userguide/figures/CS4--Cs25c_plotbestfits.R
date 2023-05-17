###################################################
### code chunk number 36: Cs25c_plotbestfits
###################################################
# plot the fits
ylbl <- rownames(dat.z)
w_ts <- seq(dim(dat.z)[2])
par(mfcol = c(3, 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:N.ts) {
  up <- fit.b$up[i, ]
  mn <- fit.b$ex[i, ]
  lo <- fit.b$lo[i, ]
  plot(w_ts, mn,
    xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", cex.lab = 1.2,
    ylim = c(min(lo), max(up))
  )
  axis(1, 12 * (0:dim(dat.z)[2]) + 1, 1980 + 0:dim(dat.z)[2])
  points(w_ts, dat.z[i, ], pch = 16, col = "blue")
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}


