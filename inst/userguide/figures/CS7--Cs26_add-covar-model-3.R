###################################################
### code chunk number 29: Cs26_add-covar-model-3
###################################################
plank.model.4 <- plank.model.3
plank.model.4$C <- matrix(list("C11", "C21", 0, 0), 4, 1)
plank.model.4$c <- d.phos


