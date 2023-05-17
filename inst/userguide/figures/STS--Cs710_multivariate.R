###################################################
### code chunk number 38: Cs710_multivariate
###################################################
dvals <- data.frame(
  x = paste0("y", 1:n),
  val = c(D, coef(fitmc, type = "matrix")$D),
  name = rep(c("true", "estimate"), each = n)
)
ggplot(dvals, aes(x = x, y = val, col = name)) +
  geom_point() +
  xlab("observation series") +
  ylab("D estimate") +
  ggtitle("D true and estimated values")


