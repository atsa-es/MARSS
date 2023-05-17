###################################################
### code chunk number 2: noshowlegend
###################################################
d <- harborSealWA
legendnames <- (unlist(dimnames(d)[2]))[2:ncol(d)]
for (i in 1:length(legendnames)) cat(paste(i, legendnames[i], "\n", sep = " "))


