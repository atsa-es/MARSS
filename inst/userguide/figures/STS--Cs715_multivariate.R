###################################################
### code chunk number 43: Cs715_multivariate
###################################################
statesm2 <- tsSmooth(fitm2)
statesm2$name <- "multiple w different Rs"
df <- rbind(true, statesu, statesm, statesmc, statesm2)

ggplot(df, aes(x = t, y = .estimate, col = name)) +
  geom_line() +
  facet_wrap(~.rownames, scale = "free_y")


