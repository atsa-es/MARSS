###################################################
### code chunk number 4: Cs002_fig1
###################################################
# Code for plotting raw Okanagan redd counts
plot(okanaganRedds[, 1], okanaganRedds[, 2],
  xlab = "", ylab = "Redd counts", main = "", col = "red", pch = 1
)
points(okanaganRedds[, 1], okanaganRedds[, 3], col = "blue", pch = 2)
legend("topleft",
  inset = 0.1, legend = c("Aerial survey", "Ground survey"),
  col = c("red", "blue"), pch = c(1, 2)
)


