###################################################
### code chunk number 23: Cs503_fitted
###################################################
ggplot(fitted2, aes(x = t, y = .estimate)) +
  geom_line() +
  facet_wrap(~.rownames, ncol = 1, scale = "free_y")


