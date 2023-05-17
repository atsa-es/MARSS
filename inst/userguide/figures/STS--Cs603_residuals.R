###################################################
### code chunk number 28: Cs603_residuals
###################################################
df2 <- data.frame(t = as.numeric(time(resids1)), resids = resids2$.std.resids, name = "MARSS")
df1 <- data.frame(t = df2$t, resids = as.numeric(resids1), name = "StructTS")
df3 <- data.frame(t = df2$t, resids = df1$resids - df2$resids, name = "difference")
df <- rbind(df1, df2, df3)

ggplot(df, aes(x = t, y = resids)) +
  geom_line() +
  facet_wrap(~name, ncol = 1, scale = "free_y") +
  ggtitle("Marginal standardized model residuals")


