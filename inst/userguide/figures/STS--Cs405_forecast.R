###################################################
### code chunk number 20: Cs405_forecast
###################################################
fr1 <- forecast:::forecast.StructTS(fit1, h = 10)
fr2 <- forecast(fit2, h = 10)
p1 <- ggplot2::autoplot(fr1, include = 8)
p2 <- ggplot2::autoplot(fr2, include = 8)
gridExtra::grid.arrange(p1, p2, nrow = 1)


