###################################################
### code chunk number 3: fig1
###################################################
d <- harborSealWA
dat <- d[, 2:ncol(d)] # first col is years
x <- d[, 1] # first col is years
n <- ncol(dat) # num time series

# set up the graphical parameters to give each data a unique line, color and width
options(warn = -99)
ltys <- matrix(1, nrow = n)
cols <- matrix(1:4, nrow = n)
lwds <- matrix(1:2, nrow = n)
pchs <- matrix(as.character(c(1:n)), nrow = n)
options(warn = 0)

graphics::matplot(x, dat, xlab = "", ylab = "log(counts)", type = "b", pch = pchs, lty = ltys, col = cols, lwd = lwds, bty = "L")
title("Puget Sound Harbor Seal Surveys")


