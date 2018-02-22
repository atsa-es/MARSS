# ###########################################
# This compares output from two different MARSS versions using the R code in the doc folder
# How to run
# Install one version of MARSS into the base R library R_HOME
# Install a second version into the local R library R_LIBS_USER
# RStudio will use R_LIBS_USER if it exists.  It does not by default so
#  you might have to create this folder to hav a local library.
# Open the unit test.R file
# RShowDoc("versiontest.R", package="MARSS")
# Change working directory to a directory where many test files can be stored (sandbox)
# Source the code.
# IMPORTANT: Using 'build and reload' from RStudio builds the package into the local
# library but does not install the doc folder (which is needed for this test)
# Use Install from zip and install from a .tar.gz file instead
# cd to folder with the source code (e.g. GitHub)
# R CMD build MARSS (build the tar.gz file)

# ###########################################

setwd("C:/Users/Eli.Holmes/Dropbox/MARSS unit tests 2018")

#make sure MARSS isn't loaded
try(detach("package:MARSS", unload=TRUE),silent=TRUE)

#One version should be in the local library
#if building from RStudio, you can set to build to local
lib.loc = Sys.getenv("R_LIBS_USER")
unittestvrs=packageVersion("MARSS", lib.loc = lib.loc)
unittestvrs #this should be new version
library(MARSS, lib.loc = lib.loc)
zscore.fun = zscore #3.9 does not have this

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
  ls.not.funs = ls()[ls()!="funs"]
  test.these = ls.not.funs[!(ls.not.funs%in%c("unittestfile","unittestfiles","unittestvrs")) & !funs]
  testNew = mget(test.these)
  save(testNew,file=paste(tag,unittestvrs,".Rdata",sep=""))
}
#detach the version
detach("package:MARSS", unload=TRUE)

#Old version of MARSS is in the R library (no local library)
lib.loc = paste(Sys.getenv("R_HOME"),"/library",sep="")
unittestvrs=packageVersion("MARSS", lib.loc = lib.loc)
unittestvrs
library(MARSS, lib.loc = lib.loc)
cat("\n\nRunning code with MARSS version", as.character(unittestvrs), "\n")
for(unittestfile in unittestfiles){
  rm(list = ls()[!(ls()%in%c("unittestfile","unittestfiles","unittestvrs","zscore.fun"))])
  tag=strsplit(unittestfile,"/")[[1]]
  tag=tag[length(tag)]
  tag=strsplit(tag,"[.]")[[1]][1]
  if(!exists(zscore)){zscore=zscore.fun}
  cat("Running ",unittestfile, "\n")
  sink(paste("outputOld-",tag,".txt",sep=""))
  set.seed(10)
  try(source(unittestfile))
  sink()
  funs=sapply(ls(),function(x){isTRUE(class(get(x))=="function")})
  ls.not.funs = ls()[ls()!="funs"]
  test.these = ls.not.funs[!(ls.not.funs%in%c("unittestfile","unittestfiles","unittestvrs")) & !funs]
  testOld = mget(test.these)
  save(testOld,file=paste(tag,unittestvrs,".Rdata",sep=""))
}
detach("package:MARSS", unload=TRUE)

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