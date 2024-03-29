###################################################
### code chunk number 16: Cs14_disfig2
###################################################
# Compare to the distance traveled per day if you used the raw data
distance.noerr <- array(NA, dim = c(dim(dat)[2] - 1, 1))
for (i in 2:dim(dat)[2]) {
  distance.noerr[i - 1] <- GCDF(dat[1, i - 1], dat[1, i], dat[2, i - 1], dat[2, i])
}
hist(distance.noerr) # make a histogram of distance traveled per day


