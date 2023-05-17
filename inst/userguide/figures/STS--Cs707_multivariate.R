###################################################
### code chunk number 35: Cs707_multivariate
###################################################
true <- data.frame(
  .rownames = rep(c("X1", "X2"), each = TT), t = t,
  .estimate = c(level, trend), .se = NA, name = "true"
)
statesu <- tsSmooth(fitu)
statesu$name <- "one bad"
statesm <- tsSmooth(fitm)
statesm$name <- "multiple bad"
df <- rbind(true, statesu, statesm)

ggplot(df, aes(x = t, y = .estimate, col = name)) +
  geom_line() +
  facet_wrap(~.rownames, scale = "free_y")


