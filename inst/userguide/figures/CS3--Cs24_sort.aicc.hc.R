###################################################
### code chunk number 34: Cs24_sort.aicc.hc
###################################################
min.AICc <- order(out.tab.hc$AICc)
out.tab.hc <- out.tab.hc[min.AICc, ]
out.tab.hc$delta.AICc <- out.tab.hc$AICc - out.tab.hc$AICc[1]
out.tab.hc$rel.like <- exp(-1 * out.tab.hc$delta.AICc / 2)
out.tab.hc$AIC.weight <- out.tab.hc$rel.like / sum(out.tab.hc$rel.like)


