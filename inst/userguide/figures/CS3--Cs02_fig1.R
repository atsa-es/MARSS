###################################################
### code chunk number 10: Cs02_fig1
###################################################
par(mfrow = c(4, 3), mar = c(2, 2, 2, 2))
for (i in 2:dim(harborSeal)[2]) {
  plot(years, harborSeal[, i], xlab = "", ylab = "", 
       main = colnames(harborSeal)[i])
}


