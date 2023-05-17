###################################################
### code chunk number 20: Cs17_print-B-Ives
###################################################
# Cleaning up the B matrix for printing
B.Ives.ML <- matrix(c(.5, NA, NA, NA, -.39, .076, NA, .1, NA, -.02, .77, NA, NA, -.1, NA, .55), 4, 4)
B.Ives.Obs <- matrix(c(.48, NA, NA, NA, -.39, .25, NA, .1, NA, -.17, .74, 0, NA, -.11, 0, .6), 4, 4)
B.Ives <- B.Ives.Obs
rownames(B.Ives) <- colnames(B.Ives) <- c("LP", "SP", "D", "ND")
print(B.Ives, digits = 2, na.print = "--")


