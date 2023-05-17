###################################################
### code chunk number 63: Cs705_globaltemp
###################################################
out_temp <- KFS(kfas_temp_stoch$model)
df <- data.frame(
  t = as.vector(time(coef(out_temp))),
  KFAS = coef(out_temp),
  `MARSS BFGS` = tsSmooth(marss_temp_stoch_bfgs, type = "xtT")$.estimate,
  `MARSS EM` = tsSmooth(marss_temp_stoch_bfgs, type = "xtT")$.estimate
)
df <- pivot_longer(df, c(KFAS, MARSS.BFGS, MARSS.EM), names_to = "model", values_to = "value")
ggplot(df, aes(x = t, y = value)) +
  geom_line() +
  facet_wrap(~model)


