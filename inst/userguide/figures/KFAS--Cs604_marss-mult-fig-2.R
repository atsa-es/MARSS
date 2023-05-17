###################################################
### code chunk number 56: Cs604_marss-mult-fig-2
###################################################
require(ggplot2)
df1 <- as.data.frame(conf_kfas_NA)
df1$name <- "KFAS"
df2 <- conf_marss_NA[, c("estimate", "Lo 95", "Hi 95")]
df2$name <- "MARSS"
colnames(df2) <- colnames(df1)
df <- rbind(df1, df2)
df$t <- as.vector(time(NileNA))
df$y <- conf_marss_NA$y
ggplot(df, aes(x = t, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey") +
  geom_line() +
  ylab("Predicted Annual Flow") +
  xlab("") +
  ggtitle("River Nile with 95% CIs on estimate") +
  facet_wrap(~name)


