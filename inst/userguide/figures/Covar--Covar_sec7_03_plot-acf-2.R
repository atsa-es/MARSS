###################################################
### code chunk number 26: Covar_sec7_03_plot-acf-2
###################################################
par(mfrow = c(5, 2), mai = c(0.1, 0.5, 0.2, 0.1), omi = c(0.5, 0, 0, 0))
for (i in 1:5) {
  plot.ts(MARSSresiduals(seas.mod.2, type = "tt1")$model.residuals[i, ],
    ylab = "Residual", main = phytos[i], xlab = "", xaxt = "n"
  )
  abline(h = 0, lty = "dashed")
  if (i == 5) {
    axis(1,
      at = 1 + seq(0, TT - period, by = 12),
      labels = seq(fulldat[years, "Year"][1], fulldat[years, "Year"][TT])
    )
    mtext(side = 1, line = 2.7, "Time")
  }
  acf(MARSSresiduals(seas.mod.2, type = "tt1")$model.residuals[i, ], lag.max = period, na.action = na.pass)
  if (i == 5) {
    axis(1, at = c(0, seq(period)))
    mtext(side = 1, line = 2.7, "Time lag")
  }
}


