###################################################
### code chunk number 19: Cs_15_plot-dlm-ACF
###################################################
# use layout to get nicer plots
layout(matrix(c(0, 1, 1, 1, 0), 1, 5, byrow = TRUE))
# set up plotting space
par(mar = c(4, 4, 1, 0), oma = c(0, 0, 0, 0.5))
# ACF of innovations
acf(t(innov), lwd = 2, lag.max = 10)


