###################################################
### code chunk number 12: Cs10_figwolfcov
###################################################
cor.fun <- function(x, y) {
  text(0.5, 0.5, format(cor(x, y), digits = 2), cex = 2)
}
pairs(t(z.score.clim.dat), lower.panel = cor.fun)


