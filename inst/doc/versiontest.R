# ###########################################
# This compares output from two different MARSS versions using the R code in the doc folder
# How to run
# Install one version of MARSS into the base R library
# Install a second version into the local R library
# Open the unit test.R file
# RShowDoc("versiontest.R", package="MARSS")
# Change working directory to a directory where many test files can be stored (sandbox)
# Source the code.
# Note: Using 'build and reload' from RStudio builds the package into the local
# library but does not install the doc or help files
# Use Install from zip and install from a .tar.gz file instead
# ###########################################

#make sure MARSS isn't loaded
try(detach(package:MARSS),silent=TRUE)

#New version should be in the local library
lib.loc = Sys.getenv("R_LIBS_USER")
unittestvrs=packageVersion("MARSS", lib.loc = lib.loc)
library(MARSS, lib.loc = lib.loc)

#Get whatever code files are in the doc directory; these are tested
unittestfiles = dir(path=paste(lib.loc,"/MARSS/doc",sep=""), pattern="*[.]R$", full.names = TRUE)
unittestfiles = unittestfiles[unittestfiles!=paste(lib.loc,"/MARSS/doc/versiontest.R",sep="")]

cat("Running code with MARSS version", as.character(unittestvrs), "\n")
for(unittestfile in unittestfiles){
  #clean the workspace but keep objects needed for the unit test
  rm(list = ls()[!(ls()%in%c("unittestfile","unittestfiles","unittestvrs"))])
  #set up name for log files
  tag=strsplit(unittestfile,"/")[[1]]
  tag=tag[length(tag)]
  tag=strsplit(tag,"[.]")[[1]][1]
  #run the code which will create objects
  cat("Running ",unittestfile, "\n")
  sink(paste("outputNew-",tag,".txt",sep=""))
  #wrapped in try so it keeps going if the code has a problem
  #set the seed so any random nums are the same
  set.seed(10)
  try(source(unittestfile))
  sink()
  #make a list of objects created by the test code
  funs=sapply(ls(),function(x){isTRUE(class(get(x))=="function")})
  test.these = ls()[!(ls()%in%c("unittestfile","unittestfiles","unittestvrs")) & !funs]
  testNew = mget(test.these)
  save(testNew,file=paste(tag,unittestvrs,".Rdata",sep=""))
}
#detach the new version
detach(package:MARSS)

#Repeat for an older version of MARSS which is in the R library (no local library)
lib.loc = paste(Sys.getenv("R_HOME"),"/library",sep="")
unittestvrs=packageVersion("MARSS", lib.loc = lib.loc)
library(MARSS, lib.loc = lib.loc)
cat("\n\nRunning code with MARSS version", as.character(unittestvrs), "\n")
for(unittestfile in unittestfiles){
  rm(list = ls()[!(ls()%in%c("unittestfile","unittestfiles","unittestvrs"))])
  tag=strsplit(unittestfile,"/")[[1]]
  tag=tag[length(tag)]
  tag=strsplit(tag,"[.]")[[1]][1]
  cat("Running ",unittestfile, "\n")
  sink(paste("outputOld-",tag,".txt",sep=""))
  set.seed(10)
  try(source(unittestfile))
  sink()
  funs=sapply(ls(),function(x){isTRUE(class(get(x))=="function")})
  test.these = ls()[!(ls()%in%c("unittestfile","unittestfiles","unittestvrs")) & !funs]
  testOld = mget(test.these)
  save(testOld,file=paste(tag,unittestvrs,".Rdata",sep=""))
}
detach(package:MARSS)

#Now start comparing the lists made using different versions of MARSS
cat("\n\nStarting object comparisons\n")
for(unittestfile in unittestfiles){
  #Get the file name
  tag=strsplit(unittestfile,"/")[[1]]
  tag=tag[length(tag)]
  tag=strsplit(tag,"[.]")[[1]][1]
  #Load in the 2 lists, testNew and testOld
  vrs=packageVersion("MARSS", lib.loc = Sys.getenv("R_LIBS_USER"))
  load(file=paste(tag,vrs,".Rdata",sep=""))
  lib.loc = paste(Sys.getenv("R_HOME"),"/library",sep="")
  vrs=packageVersion("MARSS", lib.loc = lib.loc)
  load(file=paste(tag,vrs,".Rdata",sep=""))
  
  #Compare the lists and report any differences
  cat("Checking ", tag, "\n")
  if(!identical(names(testNew), names(testOld))){
    cat("ERROR: Names of the test lists not identical\n\n")
    next
  }
  good=rep(TRUE,length(names(testNew)))
  for(ii in 1:length(names(testNew))){
    if(!identical(testNew[[ii]], testOld[[ii]])) good[ii] = FALSE
  }
  if(!all(good)){
    cat("ERROR: The following objects are not identical\n")
    cat(names(testNew)[!good])
    cat("\n\n")
  }else{
    cat("PASSED\n\n")
  }
}