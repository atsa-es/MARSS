###################################################
### code chunk number 7: Cs05_badtagfig
###################################################
# load the map package; you have to install it first
library(maps)
# Read in our noisy data (no missing values)
pdat <- loggerheadNoisy # for plotting
turtlename <- "BigMama"
par(mai = c(0, 0, 0, 0), mfrow = c(1, 1))
map("state", region = c(
  "florida", "georgia", "south carolina", "north carolina",
  "virginia", "delaware", "new jersey", "maryland"
), xlim = c(-85, -70))
points(pdat$lon[which(pdat$turtle == turtlename)], pdat$lat[which(pdat$turtle == turtlename)],
  col = "blue", pch = 21, cex = 0.7
)


