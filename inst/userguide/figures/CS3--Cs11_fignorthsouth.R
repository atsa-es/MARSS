###################################################
### code chunk number 20: Cs11_fignorthsouth
###################################################
best.fit <- fits[min.AICc][[1]]
graphics::matplot(years, t(best.fit$states - best.fit$states[, 1]),
  ylab = "Abundance index", xlab = "",
  type = "l", lwd = 2, col = "black"
)
legend("topleft", c("North Coastal", "Inland Straits", "Puget Sound", "South Coastal"), lwd = 2, lty = c(1:4), bty = "n")


