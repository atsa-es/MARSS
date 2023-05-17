###################################################
### code chunk number 6: Cs104_structTS-level
###################################################
require(tidyr)
require(ggplot2)
df1 <- as.data.frame(fit1$fitted)
vars <- colnames(df1)
df2 <- as.data.frame(t(fit3$kf$xtt))
colnames(df2) <- vars
df3 <- as.data.frame(t(fit4$kf$xtt))
colnames(df3) <- vars
df1$model <- "StructTS"
df2$model <- "MARSS BFGS"
df3$model <- "MARSS EM"
df1$t <- as.vector(time(fit1$fitted))
df2$t <- df1$t
df3$t <- df1$t
df <- rbind(df1, df2, df3)
df <- df %>% pivot_longer(all_of(vars))
ggplot(df, aes(x = t, y = value)) +
  geom_line() +
  facet_wrap(~model) +
  ggtitle("Level estimate from model fit with StructTS and MARSS")


