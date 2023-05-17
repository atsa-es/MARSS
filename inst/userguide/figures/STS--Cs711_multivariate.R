###################################################
### code chunk number 39: Cs711_multivariate
###################################################
statesmc <- tsSmooth(fitmc)
statesmc$name <- "multiple w covariate"
df <- rbind(true, statesu, statesm, statesmc)

ggplot(df, aes(x = t, y = .estimate, col = name)) +
  geom_line() +
  facet_wrap(~.rownames, scale = "free_y")


