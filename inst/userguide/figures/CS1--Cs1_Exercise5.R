###################################################
### code chunk number 27: Cs1_Exercise5
###################################################
# If you have your data in a tab delimited file with a header
# This is how you would read it in using file.choose()
# to call up a directory browser.
# However, the package has the datasets for the examples
# dat=read.table(file.choose(), skip=1)
# dat=as.matrix(dat)
dat <- wilddogs
CSEGriskfigure(dat, CI.method = "hessian", silent = TRUE)


