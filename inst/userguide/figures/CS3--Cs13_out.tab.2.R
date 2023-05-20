###################################################
### code chunk number 23: Cs13_out.tab.2
###################################################
min.AICc <- order(out.tab$AICc)
out.tab.2 <- out.tab[min.AICc, ]
fits <- fits[min.AICc]
out.tab.2$delta.AICc <- out.tab.2$AICc - out.tab.2$AICc[1]
out.tab.2$rel.like <- exp(-1 * out.tab.2$delta.AICc / 2)
out.tab.2$AIC.weight <- out.tab.2$rel.like / sum(out.tab.2$rel.like)


