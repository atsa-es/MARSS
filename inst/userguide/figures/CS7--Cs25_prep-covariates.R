###################################################
### code chunk number 28: Cs25_prep-covariates
###################################################
# transpose to make time go across columns
# drop=FALSE so that R doesn't change our matrix to a vector
phos <- t(log(ivesDataByWeek[, "Phosph", drop = FALSE]))
d.phos <- (phos - apply(phos, 1, mean, na.rm = TRUE))


