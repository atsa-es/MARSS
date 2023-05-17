###################################################
### code chunk number 65: Cs707_globaltemp
###################################################
kfs <- KFS(kfas_temp_stoch$model, smoothing = "disturbance")
resid_kfas <- rstandard(kfs,
  type = "state",
  standardization_type = "cholesky"
)
resid_marss <- residuals(marss_temp_stoch_bfgs,
  type = "tT",
  standardization = "Block.Cholesky"
)
resid_marss <- subset(resid_marss, name == "state")
resid_test <- residuals(marss_test,
  type = "tT",
  standardization = "Block.Cholesky"
)
resid_test <- subset(resid_test, name == "state")

df <- data.frame(
  MARSS = resid_marss$.std.resids,
  KFAS = as.vector(resid_kfas),
  MARSS.test = resid_test$.std.resids
)
df$diff.est <- df$MARSS - df$KFAS
df$diff.id <- df$MARSS.test - df$KFAS
df$t <- as.vector(time(kfas_temp_stoch$model$y))
df$name <- "state"
df1 <- pivot_longer(df, c(MARSS, MARSS.test, KFAS, diff.est, diff.id), names_to = "model", values_to = "value")

kfs <- KFS(kfas_temp_stoch$model)
resid_kfas <- residuals(kfs, type = "pearson")
resid_marss <- MARSSresiduals(marss_temp_stoch_bfgs, type = "tT")
resid_test <- MARSSresiduals(marss_test, type = "tT")
df <- rbind(
  cbind(as.data.frame(t(resid_marss$model.residuals)), model = "MARSS") %>% pivot_longer(c("HL", "Folland")),
  cbind(as.data.frame(t(resid_test$model.residuals)), model = "MARSS.test") %>% pivot_longer(c("HL", "Folland")),
  cbind(as.data.frame(resid_kfas), model = "KFAS") %>% pivot_longer(c("HL", "Folland")),
  cbind(as.data.frame(t(resid_marss$model.residuals) - resid_kfas), model = "diff.est") %>% pivot_longer(c("HL", "Folland")),
  cbind(as.data.frame(t(resid_test$model.residuals) - resid_kfas), model = "diff.id") %>% pivot_longer(c("HL", "Folland"))
)
df$t <- rep(as.vector(time(kfas_temp_stoch$model$y)), 2 * 5)

df <- rbind(df, df1[, colnames(df)])

ggplot(subset(df, model %in% c("diff.est", "diff.id")), aes(x = t, y = value, col = model)) +
  geom_line(na.rm = TRUE) +
  facet_wrap(~name) +
  xlab("") +
  ggtitle("Difference in residuals KFAS vs MARSS")


