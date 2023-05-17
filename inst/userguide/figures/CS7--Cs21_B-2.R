###################################################
### code chunk number 24: Cs21_B-2
###################################################
B.2 <- matrix(list(0), 4, 4) # set up the list matrix
diag(B.2) <- c("B11", "B22", "B33", "B44") # give names to diagonals
# and names to the estimated non-diagonals
B.2[1, 2] <- "B12"
B.2[2, 3] <- "B23"
B.2[2, 4] <- "B24"
B.2[4, 2] <- "B42"
print(B.2)


