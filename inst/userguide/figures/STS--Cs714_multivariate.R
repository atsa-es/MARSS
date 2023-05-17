###################################################
### code chunk number 42: Cs714_multivariate
###################################################
rvals <- data.frame(
  x = paste0("y", 1:n),
  val = c(r2, coef(fitm2)$R),
  name = rep(c("true", "estimate"), each = n)
)
ggplot(rvals, aes(x = x, y = val, col = name)) +
  geom_point() +
  xlab("observation series") +
  ylab("R variance estimate") +
  ggtitle("R true and estimated values")


