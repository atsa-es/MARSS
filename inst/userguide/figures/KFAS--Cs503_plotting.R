###################################################
### code chunk number 51: Cs503_plotting
###################################################
require(ggplot2)
df <- cbind(conf_marss1, pred_marss1[, c(".lwr", ".upr")])
ggplot(df, aes(x = t, y = .fitted)) +
  geom_ribbon(aes(ymin = .lwr, ymax = .upr), fill = "grey") +
  geom_ribbon(aes(ymin = .conf.low, ymax = .conf.up), fill = "blue", alpha = 0.25) +
  geom_line(linetype = 2) +
  ylab("Predicted Annual Flow") +
  xlab("") +
  ggtitle("River Nile")


