###################################################
### code chunk number 58: Cs606_missing-values
###################################################
require(ggplot2)
require(tidyr)
df <- rbind(fitted_kfas_NA, fitted_marss_NA)
df$t <- as.vector(time(NileNA))
df$y <- conf_marss_NA$y
df <- tidyr::pivot_longer(df, c(smooth, one.step.ahead), names_to = "type", values_to = "value")
ggplot(df, aes(x = t, y = value, col = type)) +
  geom_point(aes(x = t, y = y), col = "blue", size = 0.5, na.rm = TRUE) +
  geom_line() +
  ylab("Predicted Annual Flow") +
  xlab("") +
  ggtitle("River Nile - smoothed versus filtered") +
  facet_wrap(~name, ncol = 1)


