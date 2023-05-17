###################################################
### code chunk number 35: Cs25_out.tab.2
###################################################
out.tab.hc$AIC.weight <- round(out.tab.hc$AIC.weight, digits = 3)
out.tab.hc$delta.AICc <- round(out.tab.hc$delta.AICc, digits = 2)
print(out.tab.hc[, c("H", "Q", "delta.AICc", "AIC.weight")], row.names = FALSE)


