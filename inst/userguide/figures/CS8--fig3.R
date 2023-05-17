###################################################
### code chunk number 18: fig3
###################################################
d <- rockfish
dat <- d[, 2:ncol(d)] # first col is years
x <- d[, 1] # first col is years
n <- nrow(dat) # num time series

# set up the graphical parameters to give each data a unique line, color and width
options(warn = -99)
ltys <- matrix(1, nrow = n)
cols <- matrix(1:4, nrow = n)
lwds <- matrix(1:2, nrow = n)
pchs <- matrix(as.character(c(1:n)), nrow = n)
options(warn = 0)

meany <- matrix(apply(dat, 2, mean, na.rm = T), nrow = nrow(dat), ncol = ncol(dat), byrow = T) # take off the mean; for better plotting and to mask pattern
adj <- matrix(c(1.5, 0, -2, 1, 2, -2, -1, 0, 2), nrow = nrow(dat), ncol = ncol(dat), byrow = T) # just to mask pattern in data
graphics::matplot(x, dat - meany + adj, xlab = "", ylab = "Index of log(cpue)", type = "b", pch = pchs, lty = ltys, col = cols, lwd = lwds)
title("Puget Sound Total Rockfish Indices")


