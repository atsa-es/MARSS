###################################################
### code chunk number 33: Cs24_plottrends
###################################################
# get ts of trends
ts.trends <- t(trends.rot)
par(mfrow = c(ceiling(dim(ts.trends)[2] / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
# loop over each trend
for (i in 1:dim(ts.trends)[2]) {
  # set up plot area
  plot(ts.trends[, i],
    ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
    type = "n", lwd = 2, bty = "L",
    xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )
  # draw zero-line
  abline(h = 0, col = "gray")
  # plot trend line
  par(new = TRUE)
  plot(ts.trends[, i],
    ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
    type = "l", lwd = 2, bty = "L",
    xlab = "", ylab = "", xaxt = "n"
  )
  # add panel labels
  mtext(paste("Trend", i, sep = " "), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(dat.spp.1980)[2]) + 1, 1980 + 0:dim(dat.spp.1980)[2])
} # end i loop (trends)


