###################################################
### code chunk number 24: Cs14_out.tab.2
###################################################
out.tab.2$AIC.weight <- round(out.tab.2$AIC.weight, digits = 3)
out.tab.2$delta.AICc <- round(out.tab.2$delta.AICc, digits = 2)
print(out.tab.2[1:10, c("H", "Q", "delta.AICc", "AIC.weight")], row.names = FALSE)


