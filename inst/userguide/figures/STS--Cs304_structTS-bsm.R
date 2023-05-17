###################################################
### code chunk number 15: Cs304_structTS-bsm
###################################################
require(tidyr)
require(ggplot2)
df1 <- as.data.frame(fit1$fitted)
vars <- colnames(df1)
df2 <- as.data.frame(t(fit3$kf$xtt)[, 1:3])
colnames(df2) <- vars
df3 <- as.data.frame(t(fit4$kf$xtt)[, 1:3])
colnames(df3) <- vars
df1$model <- "StructTS"
df2$model <- "MARSS BFGS"
df3$model <- "MARSS EM"
df1$t <- as.vector(time(fit1$fitted))
df1$Qtr <- as.vector(cycle(fit1$fitted))
df2$t <- df1$t
df2$Qtr <- df1$Qtr
df3$t <- df1$t
df3$Qtr <- df1$Qtr
df <- rbind(df1, df2, df3)
df <- subset(df, Qtr == 1) %>% pivot_longer(all_of(vars))
ggplot(df, aes(x = t, y = value, color = model, linetype = model, shape = model)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free") +
  scale_linetype_manual("model", values = c(1, 1, 0)) +
  scale_shape_manual("model", values = c(NA, NA, 15))


