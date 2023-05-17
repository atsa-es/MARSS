###################################################
### code chunk number 16: Cs_13_plot-dlmQQ
###################################################
# use layout to get nicer plots
layout(matrix(c(0, 1, 1, 1, 0), 1, 5, byrow = TRUE))
# set up L plotting space
par(mar = c(4, 4, 1, 0), oma = c(0, 0, 0, 0.5))
# Q-Q plot of innovations
qqnorm(t(innov), main = "", pch = 16, col = "blue")
qqline(t(innov))
# set up R plotting space
# par(mar=c(4,0,1,1)) #, oma=c(0,0,0,0.5))
# boxplot of innovations
# boxplot(t(innov), axes=FALSE)


