###################################################
### code chunk number 14: Cs11_distance2
###################################################
distance <- array(NA, dim = c(dim(dat)[2] - 1, 1))
for (i in 2:dim(dat)[2]) {
  distance[i - 1] <- GCDF(
    pred.lon[i - 1], pred.lon[i],
    pred.lat[i - 1], pred.lat[i]
  )
}


