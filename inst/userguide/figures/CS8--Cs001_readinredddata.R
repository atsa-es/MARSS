###################################################
### code chunk number 3: Cs001_readinredddata
###################################################
head(okanaganRedds)
logRedds <- log(t(okanaganRedds)[c("aerial", "ground"), ])


