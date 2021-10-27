# ###########################################
# This compares output from two different MARSS versions using the R code in the doc folder
# How to run
# Install one version of MARSS into the base R library R_HOME
# Install a second version into the local R library R_LIBS_USER
# RStudio will use R_LIBS_USER if it exists.  It does not by default so
#  you might have to create this folder to have a local library.
#  look at Sys.getenv("R_LIBS_USER"). Click Install under Packages tab and see
#  where it is installing
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

if(Sys.info()['sysname']=="Windows"){
  setwd("C:/Users/Eli.Holmes/Dropbox/MARSS unit tests 2021")
  lib.new <- "C:/Program Files/R/R-3.6.2/library"
}
if(Sys.info()['sysname']=="Darwin"){
  setwd("~/Dropbox/MARSS unit tests 2021")
  lib.new <- "/Library/Frameworks/R.framework/Versions/4.1/Resources/library"
}
lib.old <- Sys.getenv("R_LIBS_USER")

# to install MARSS to correct locations if needed
# install.packages("MARSS", lib.old) #install from CRAN
# Mac: install.packages("~/Dropbox/MARSS_3.11.1.tar.gz", lib=lib.old, repos=NULL)
# Mac: install.packages("~/Dropbox/MARSS_3.11.2.tar.gz", lib=lib.new, repos=NULL)
# Win: install.packages("C:/Users/Eli.Holmes/Dropbox/MARSS_3.11.2.tar.gz", lib=lib.new, repos=NULL)

#make sure MARSS isn't loaded
try(detach("package:MARSS", unload=TRUE),silent=TRUE)

#Load new and get the R files
unittestfiles = dir(path=paste(lib.new,"/MARSS/doc",sep=""), pattern="*[.]R$", full.names = TRUE)
unittestfiles = unittestfiles[unittestfiles!=paste(path.expand(lib.new),"/MARSS/doc/versiontest.R",sep="")]

unittestvrs=packageVersion("MARSS", lib.loc = lib.new)
unittestvrs #this should be new version
library(MARSS, lib.loc = lib.new)
# #zscore.fun = zscore #3.9 does not have this
# MARSSresiduals.fun = MARSSresiduals
# MARSSresiduals_tT.fun = MARSS:::MARSSresiduals.tT
# MARSSresiduals_tt1.fun = MARSS:::MARSSresiduals.tt1
# MARSSresiduals_tt.fun = MARSS:::MARSSresiduals.tt

# file <- "AR2SS100.RData"
# if (file %in% dir("./vignettes/manual_files")) {
#   load(paste("./vignettes/manual_files/", file, sep = ""))
#   sims.exist <- TRUE
# } else {
#   sims.exist <- FALSE
# }
# file <- "CS4--model_fits.RData"
# if (file %in% dir("./vignettes/manual_files")) {
#   load(paste("./vignettes/manual_files/", file, sep = ""))
#   sims.exist <- TRUE
# } else {
#   sims.exist <- FALSE
# }

cat("Running code with MARSS version", as.character(unittestvrs), "\n")
for(unittestfile in unittestfiles){
  #clean the workspace but keep objects needed for the unit test
  rm(list = ls()[!(ls()%in%c("unittestfile","unittestfiles","unittestvrs","zscore.fun","lib.new","lib.old", "MARSSresiduals.fun", "MARSSresiduals_tT.fun", "MARSSresiduals_tt1.fun"))])
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
  funs=sapply(ls(),function(x){isTRUE(class(get(x))[1]=="function")})
  ls.not.funs = ls()[ls()!="funs"]
  test.these = ls.not.funs[!(ls.not.funs%in%c("unittestfile","unittestfiles","unittestvrs")) & !funs]
  testNew = mget(test.these)
  save(testNew,file=paste(tag,unittestvrs,".Rdata",sep=""))
}
#detach the version
detach("package:MARSS", unload=TRUE)

#Load old version of MARSS
unittestvrs=packageVersion("MARSS", lib.loc = lib.old)
unittestvrs
library(MARSS, lib.loc = lib.old)
cat("\n\nRunning code with MARSS version", as.character(unittestvrs), "\n")
for(unittestfile in unittestfiles){
  rm(list = ls()[!(ls()%in%c("unittestfile","unittestfiles","unittestvrs","zscore.fun","lib.new","lib.old", "MARSSresiduals.fun", "MARSSresiduals_tT.fun", "MARSSresiduals_tt1.fun"))])
  tag=strsplit(unittestfile,"/")[[1]]
  tag=tag[length(tag)]
  tag=strsplit(tag,"[.]")[[1]][1]
#  if(!exists("zscore")){zscore=zscore.fun}
# #  if(!exists("MARSSresiduals")){
#     MARSSresiduals = MARSSresiduals.fun
#     MARSSresiduals.tT = MARSSresiduals_tT.fun
#     MARSSresiduals.tt1 = MARSSresiduals_tt1.fun
#     MARSSresiduals.tt = MARSSresiduals_tt.fun
# #  }
  cat("Running ",unittestfile, "\n")
  sink(paste("outputOld-",tag,".txt",sep=""))
  set.seed(10)
  testfile <- unittestfile
  if(unittestvrs=="3.10.14" && tag %in% c("Chapter_DFA")){
    con <- file(unittestfile)
    sourcercode <- readLines(con)
    close(con)
    sourcercode <- str_replace(sourcercode, '"tT', '"smoothations')
    sourcercode <- str_replace(sourcercode, 'name==', 'type==')
    con <- file('testfile.R')
    writeLines(sourcercode, con)
    testfile <- 'testfile.R'
    close(con)
    rm(sourcercode, con)
  }
  if(unittestvrs=="3.10.14" && tag %in% c("Chapter_SealTrend", "Chapter_StructuralBreaks")){
    con <- file(unittestfile)
    sourcercode <- readLines(con)
    close(con)
    sourcercode <- str_replace(sourcercode, 'type = \"tt1\"', 'conditioning = \"t-1\"')
    sourcercode <- str_replace(sourcercode, 'type = \"tT\"', 'conditioning = \"T\"')
    sourcercode <- str_replace(sourcercode, '[$]model[.]residuals', '$residuals')
    con <- file('testfile.R')
    writeLines(sourcercode, con)
    testfile <- 'testfile.R'
    close(con)
    rm(sourcercode, con)
  }
  if(unittestvrs=="3.10.14" && tag=="Chapter_Structural_TS"){
    con <- file(unittestfile)
    sourcercode <- readLines(con)
    close(con)
    sourcercode <- str_replace(sourcercode, 'forecast::forecast[(]fit2', 'forecast.marssMLE(fit2')
    con <- file('testfile.R')
    writeLines(sourcercode, con)
    testfile <- 'testfile.R'
    close(con)
    rm(sourcercode, con)
  }
  if(unittestvrs=="3.10.14" && tag=="Quick_Examples"){
    con <- file(unittestfile)
    sourcercode <- readLines(con)
    close(con)
    sourcercode <- str_replace(sourcercode, 'tsSmooth[(]kemfit[)]', 'NULL')
    con <- file('testfile.R')
    writeLines(sourcercode, con)
    testfile <- 'testfile.R'
    close(con)
    rm(sourcercode, con)
  }
  try(source(testfile, echo=TRUE))
  rm(testfile)
  sink()
  funs=sapply(ls(),function(x){isTRUE(class(get(x))[1]=="function")})
  ls.not.funs = ls()[ls()!="funs"]
  test.these = ls.not.funs[!(ls.not.funs%in%c("unittestfile","unittestfiles","unittestvrs")) & !funs]
  testOld = mget(test.these)
  save(testOld,file=paste(tag,unittestvrs,".Rdata",sep=""))
}
detach("package:MARSS", unload=TRUE)

#Now start comparing the lists made using different versions of MARSS
cat("\n\nStarting object comparisons\n")
# for 3.10.14 since many attribute diffs
idfun <- function(x, y){ isTRUE(all.equal(x, y, check.attributes = FALSE)) }
# for 3.11.1
#idfun <- function(x, y){ identical(x, y) }
for(unittestfile in unittestfiles){
  #Get the file name
  tag=strsplit(unittestfile,"/")[[1]]
  tag=tag[length(tag)]
  tag=strsplit(tag,"[.]")[[1]][1]
  #Load in the 2 lists, testNew and testOld
  vrs=packageVersion("MARSS", lib.loc = lib.new)
  load(file=paste(tag,vrs,".Rdata",sep=""))
  vrs=packageVersion("MARSS", lib.loc = lib.old)
  load(file=paste(tag,vrs,".Rdata",sep=""))
  
  #Compare the lists and report any differences
  cat("Checking ", tag, "\n")
  if(!idfun(names(testNew), names(testOld))){
    cat("ERROR: Names of the test lists not identical\n")
    cat("testNew has these not in testOld", setdiff(names(testNew), names(testOld)), "\n")
    cat("testOld has these not in testNew", setdiff(names(testOld), names(testNew)), "\n\n")
    next
  }
  good=rep(TRUE,length(names(testNew)))
  for(ii in 1:length(names(testNew))){
    if( names(testNew)[ii] %in% c("bfgsfit.time", "kemfit.time")) next
    if(inherits(testNew[[ii]], "marssMLE")){
      attr(testNew[[ii]][["call"]][["data"]], "model.tsp") <- NULL
      attr(testNew[[ii]][["model"]], "model.tsp") <- NULL
      attr(testNew[[ii]][["marss"]], "model.tsp") <- NULL
      for(kk in c("model", "marss")){
        attr(testNew[[ii]][[kk]], "equation") <- NULL
        attr(testOld[[ii]][[kk]], "equation") <- NULL
      }
      if(inherits(testNew[[ii]]$call$inits, "marssMLE")){
        attr(testNew[[ii]][["call"]][["inits"]][["call"]][["data"]], "model.tsp") <- NULL
        for(kk in c("model", "marss")){
          attr(testNew[[ii]]$call$inits[[kk]], "equation") <- NULL
          attr(testNew[[ii]]$call$inits[[kk]], "model.tsp") <- NULL
          attr(testOld[[ii]]$call$inits[[kk]], "equation") <- NULL
        }
      }
    }
    if(inherits(testNew[[ii]], "marssMODEL")){
      attr(testNew[[ii]], "equation") <- NULL
      attr(testNew[[ii]], "model.tsp") <- NULL
      attr(testOld[[ii]], "equation") <- NULL
    }
    if(inherits(testNew[[ii]], "SSModel")){
      attr(testNew[[ii]]$terms, ".Environment") <- NULL
      attr(testOld[[ii]]$terms, ".Environment") <- NULL
    }
    if(inherits(testNew[[ii]], "gg")){
      next
    }
    if(inherits(testNew[[ii]], "list")){
      for(iii in 1:length(testNew[[ii]])){
        if(inherits(testNew[[ii]][[iii]], "SSModel")){
            attr(testNew[[ii]][[iii]]$terms, ".Environment") <- NULL
            attr(testOld[[ii]][[iii]]$terms, ".Environment") <- NULL
          }
          if(inherits(testNew[[ii]][[iii]], "marssMLE")){
            attr(testNew[[ii]][[iii]][["call"]][["data"]], "model.tsp") <- NULL
            for(kk in c("model", "marss")){
            attr(testNew[[ii]][[iii]][[kk]], "equation") <- NULL
            attr(testNew[[ii]][[iii]][[kk]], "model.tsp") <- NULL
            attr(testOld[[ii]][[iii]][[kk]], "equation") <- NULL
          }
          if(inherits(testNew[[ii]][[iii]]$call$inits, "marssMLE")){
            attr(testNew[[ii]][[iii]][["call"]][["inits"]][["call"]][["data"]], "model.tsp") <- NULL
            for(kk in c("model", "marss")){
              attr(testNew[[ii]][[iii]]$call$inits[[kk]], "equation") <- NULL
              attr(testNew[[ii]][[iii]]$call$inits[[kk]], "model.tsp") <- NULL
              attr(testOld[[ii]][[iii]]$call$inits[[kk]], "equation") <- NULL
            }
          }
        }
        if(inherits(testNew[[ii]][[iii]], "marssMODEL")){
          attr(testNew[[ii]][[iii]], "equation") <- NULL
          attr(testNew[[ii]][[iii]], "model.tsp") <- NULL
          attr(testOld[[ii]][[iii]], "equation") <- NULL
        }
      }
    }
    if(!idfun(testNew[[ii]], testOld[[ii]])){
      good[ii] = FALSE
      if(inherits(testNew[[ii]], "marssMLE")){
        for(iii in names(testNew[[ii]][["par"]])){
          if(iii %in% c("G","H","L")) next
          if(!idfun(testNew[[ii]][["par"]][iii], testOld[[ii]][["par"]][iii])){
            cat("Warning:", names(testNew)[ii],"par",iii,"not identical\n")
          }else{
            #cat(names(testNew)[ii],"par",iii,"idfun\n")
          }
        }
        for(iii in names(testNew[[ii]][["call"]])){
          if(!idfun(testNew[[ii]][["call"]][iii], testOld[[ii]][["call"]][iii])){
            cat("Warning:", names(testNew)[ii],"call",iii,"not identical\n")
          }else{
            #cat(names(testNew)[ii],"call",iii,"idfun\n")
          }
        }
        for(iii in names(testNew[[ii]])){
          if (!(iii %in% names(testOld[[ii]]))) {
            cat("Warning:", names(testNew)[ii],iii,"not in testOld.\n")
          } else {
          if(!idfun(testNew[[ii]][iii], testOld[[ii]][iii])){
            cat("Warning:", names(testNew)[ii],iii,"not identical\n")
          }
            }
        }
      }
      if(inherits(testNew[[ii]], "matrix") || inherits(testNew[[ii]], "array")){
        if(!idfun(dim(testNew[[ii]]), dim(testOld[[ii]]))){
          cat("Warning: dims of", names(testNew)[ii], "not identical\n")
          next
        }
        if(!isTRUE(all.equal(testNew[[ii]], testOld[[ii]], check.attributes = FALSE)))
          cat("Warning: values in", names(testNew)[ii], "not identical\n")
        if(!idfun(rownames(testNew[[ii]]), rownames(testOld[[ii]])))
          cat("Warning: rownames of", names(testNew)[ii], "not identical\n")
      }
      if(inherits(testNew[[ii]], "list")){
        for(kk in 1:length(testNew[[ii]])){
          if(inherits(testNew[[ii]][[kk]], "marssMLE")){
            for(iii in names(testNew[[ii]][[kk]][["par"]])){
              if(iii %in% c("G","H","L")) next
              if(!idfun(testNew[[ii]][[kk]][["par"]][iii], testOld[[ii]][[kk]][["par"]][iii])){
                cat("Warning:", names(testNew)[[ii]][kk],"par",iii,"not identical\n")
              }else{
                #cat(names(testNew)[ii],"par",iii,"idfun\n")
              }
            }
            for(iii in names(testNew[[ii]][[kk]][["call"]])){
              if(!idfun(testNew[[ii]][[kk]][["call"]][iii], testOld[[ii]][[kk]][["call"]][iii])){
                cat("Warning:", names(testNew)[[ii]][kk],"call",iii,"not identical\n")
              }else{
                #cat(names(testNew)[ii],"call",iii,"idfun\n")
              }
            }
            for(iii in names(testNew[[ii]][[kk]])){
              if(!idfun(testNew[[ii]][[kk]][iii], testOld[[ii]][[kk]][iii])){
                cat("Warning:", names(testNew)[[ii]][kk],iii,"not identical\n")
              }else{
                #cat(names(testNew)[ii],iii,"idfun\n")
              }
            }
          }
          if(inherits(testNew[[ii]][[kk]], "matrix") || inherits(testNew[[ii]][[kk]], "array")){
            if(!idfun(dim(testNew[[ii]][[kk]]), dim(testOld[[ii]][[kk]]))){
              cat("Warning: dims of", names(testNew)[ii], "$", names(testNew[[ii]])[kk], "not identical\n")
              next
            }
            if(!idfun(testNew[[ii]][[kk]], testOld[[ii]][[kk]]))
              cat("Warning: values in", names(testNew)[ii], "$", names(testNew[[ii]])[kk], "not identical\n")
            if(!idfun(rownames(testNew[[ii]][[kk]]), rownames(testOld[[ii]][[kk]])))
              cat("Warning: rownames of", names(testNew)[ii], "$", names(testNew[[ii]])[kk], "not identical\n")
            
          }
        }
      }
    }
  }
  if(!all(good)){
    cat("ERROR: The following objects are not identical\n")
    cat(names(testNew)[!good])
    cat("\n\n")
  }else{
    cat("PASSED\n\n")
  }
}