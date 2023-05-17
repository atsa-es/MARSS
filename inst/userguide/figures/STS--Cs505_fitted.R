###################################################
### code chunk number 25: Cs505_fitted
###################################################
fitted2 <- fitted(fit2, type = "ytt")
df2 <- data.frame(t = as.numeric(time(fitted1)), fitted = fitted2$.fitted, name = "MARSS")
df1 <- data.frame(t = df2$t, fitted = as.numeric(fitted1[, 1] + fitted1[, 3]), name = "StructTS")
df <- rbind(df1, df2)
df$y <- fitted2$y

ggplot(df) +
  geom_line(aes(x = t, y = fitted)) +
  geom_point(aes(x = t, y = y), col = "blue") +
  facet_wrap(~name, ncol = 1)


