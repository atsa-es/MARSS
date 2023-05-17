###################################################
### code chunk number 19: Cs10_print.table
###################################################
out.tab.1$delta.AICc <- round(out.tab.1$delta.AICc, digits = 2)
out.tab.1$AIC.weight <- round(out.tab.1$AIC.weight, digits = 3)
print(out.tab.1[, c("H", "Q", "delta.AICc", "AIC.weight")], row.names = FALSE)


