###################################################
### code chunk number 48: Cs720_multivariate
###################################################
true <- data.frame(
  .rownames = rep(c("level 1", "level 2", "trend"), each = TT), t = t,
  .estimate = c(level1, level2, trend), .se = NA, name = "true"
)
statesm1 <- tsSmooth(fitm1)
statesm1$name <- "ts 1 alone"
statesm1$.rownames[statesm1$.rownames == "X2"] <- "trend"
statesm1$.rownames[statesm1$.rownames == "X1"] <- "level 1"
statesm2 <- tsSmooth(fitm2)
statesm2$name <- "ts 2 alone"
statesm2$.rownames[statesm2$.rownames == "X2"] <- "trend"
statesm2$.rownames[statesm2$.rownames == "X1"] <- "level 2"
statesm3 <- tsSmooth(fitm3)
statesm3$name <- "ts 1 & 2 together"
statesm3$.rownames[statesm3$.rownames == "X3"] <- "trend"
statesm3$.rownames[statesm3$.rownames == "X1"] <- "level 1"
statesm3$.rownames[statesm3$.rownames == "X2"] <- "level 2"
df <- rbind(true, statesm1, statesm2, statesm3)

ggplot(df, aes(x = t, y = .estimate, col = name)) +
  geom_line() +
  facet_wrap(~.rownames, scale = "free_y", ncol = 1)


