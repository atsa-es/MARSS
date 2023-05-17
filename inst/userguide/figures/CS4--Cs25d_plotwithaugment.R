###################################################
### code chunk number 37: Cs25d_plotwithaugment
###################################################
require(ggplot2)
alpha <- 0.05
theme_set(theme_bw())
d <- residuals(the.fit, type = "tT")
d$up <- qnorm(1 - alpha / 2) * d$.sigma + d$.fitted
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted
ggplot(data = subset(d, name=="model")) +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  geom_line(aes(t, .fitted), col="blue") +
  facet_wrap(~.rownames) +
  xlab("Time Step") +
  ylab("Count")


