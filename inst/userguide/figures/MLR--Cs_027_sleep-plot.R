###################################################
### code chunk number 30: Cs_027_sleep-plot
###################################################
library(lattice)
xyplot(Reaction ~ Days | Subject, sleepstudy,
  type = c("g", "p", "r"),
  index = function(x, y) coef(lm(y ~ x))[1],
  xlab = "Days of sleep deprivation",
  ylab = "Average reaction time (ms)", aspect = "xy"
)


