# Steps for a release

Once a release is ready, all tests done and documentation checked, the Makefile in `inst/doc` is run. I usually run one target at a time starting with `make movetovignettes`. Before doing this, I build the User Guide and run it through a few compressors: Adobe and then https://github.com/atsa-es/pdfsizeopt-readme.

## Spellcheck the whole package

```
#devtools::install_github("ropensci/spelling")
spelling::spell_check_package()
allfils <- c("README.md", "DESCRIPTION", "NEWS.md")
fils <- dir("inst/userguide/manual_files", full.names = TRUE)
fils <- fils[stringr::str_detect(fils, "Rnw")]
allfils <- c(allfils, fils)
fils <- dir("inst/derivations", full.names = TRUE)
fils <- fils[stringr::str_detect(fils, "Rnw")]
allfils <- c(allfils, fils)
fils <- dir("vignettes", full.names = TRUE)
fils <- fils[stringr::str_detect(fils, "Rmd")]
allfils <- c(allfils, fils)
fils <- dir("man", full.names = TRUE)
fils <- fils[stringr::str_detect(fils, "Rd")]
allfils <- c(allfils, fils)
spelling::spell_check_files(fils, ignore = spelling::get_wordlist(pkg = "."), lang = "en_US")
```

## Derivations

Docs are in `inst/derivations`. Re-run if the derivation files change. Running the make file will update the `inst/derivations/doc` folder which is used for the files that are put in the `inst/doc` folder of a source file .tar.gz release.

```
cd inst/derivations
make
```

## User Guide

Re-run with any new release. Running the make file will update the `inst/userguide/doc` folder which is used for the files that are put in the `inst/doc` folder of a source file .tar.gz release.

```
cd inst/userguide
make
```

You can update both by 

```
cd inst
make --makefile=Makedoc
```

## Check the documentation pdfs

Will be in `inst/userguide/doc` and `inst/derivations/doc`

  - [ ] User Guide
  - [ ] EMDerivation
  - [ ] Residuals

## Compress User Guide

  - [ ] Manually compress the UserGuide.pdf in the inst/userguide/doc folder in Dropbox/MARSS/inst/doc with `https://github.com/atsa-es/pdfsizeopt-readme`
  - [ ] Rebuild package with smaller docs, and then re-Run checks on the built package to see that it passes size checks

Should be about 1.3M

## Manually update the inst/doc for GitHub

So that users install the vignettes when they install from GitHub

  - [ ] copy from Dropbox/MARSS/inst/doc into GitHub/MARSS/inst/doc
  
## Check the manual built by CRAN

```
setwd("~/Documents/GitHub")
devtools::build_manual(pkg="MARSS")
```

## Check the vignettes

Run and check and then delete the detritus.

  - [ ] Quick_Start.html
  - [ ] Learning_MARSS.html

## Basic Tests

- [ ] Run this `vignettes/versiontest.R`

- [ ] Run checks on the built package

```
rm -r ~/Dropbox/MARSS.Rcheck
R CMD check --timings --as-cran MARSS_3.11.7.tar.gz
```

## Run internal tests 

  - [ ] Source `tests/testthat/model.R` to set up models for tests

A note about these tests. I tried to break MARSS code with the tests and choose models and data that are difficult and weird. The tests will pump out many errors and warnings. Don't worry about that. Just watch if the test fails.

You can run all the tests with `devtools::test()` but I find running the tests one by one in RStudio is easier to make sense of since there are so many error messages. Open the test file and you will see a button to run tests. Some tests will show errors, but all tests should pass. Errors are correct if the test passes and result should be try-error.

These are the tests in `tests/testthat`. All should pass but some give warnings. See note `test-GDF-example` for a comment regarding differences on different OS.

- [ ] test.coef `[ FAIL 0 | WARN 8 | SKIP 0 | PASS 57 ]`

  This will warn about "MARSSparamCIs: No parSigma element returned by Hessian function.  See marssMLE object errors (MLEobj$errors)". This is ok.

- [ ] test-forecast `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 76 ]`

  There will be some errors. This is correct. Test is to make sure these error messages appear.

- [ ] test-GDF-example `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 18 ]`

  This one shows not insignificant differences in the logLik tests (fails) when run on my MacBook Air (M1, 2020) versus my old Mac (Intel) or on a Linux (Posit.Cloud R 4.3, everything updated). If run on Linux, then everything should pass as the logLik hard coded values are from that OS.


- [ ] test-KFAS `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 45 ]`

  This has the code from the `inst/manual_files/KFAS.Rnw` however the `se.fit` for output from `predict()` for KFAS models is removed because this is expected to be different. See the KFAS chapter at the end of the User Guide.


- [ ] test-MARSSkfss `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 68 ]`

  This is testing that kfss and kfas are returning the same values across a variety of models.


- [ ] test-plotting `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 160 ]`

  Lots of residuals warnings on this one as these are difficult models. "MARSSresiduals.tt1 reported warnings. See msg element of returned residuals object."" Warnings are correct, but all tests should pass. I kind of went overkill on the tests here.


- [ ] test-print `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 120 ]`

- [ ] test-residualsHarvey `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 84 ]`

  There will be some warnings about MARSSresiduals.tT reporing warnings.

- [ ] test-structTS `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 5 ]`

- [ ] test-tt `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 4403 ]`

  A long kind of excessive test but in 3.11.4 I was having trouble with the tt residuals code.

## Check any reverse depends

* [ ] `revdepcheck::revdep_check(num_workers = 4)`

## Ready for final build

```
cd inst
make
```
However I usually run it one target at a time.


Now using the final tar.gz file do the final tests.
  
## Check timings
  
```
cd ~/Dropbox
rm -r ~/Dropbox/MARSS.Rcheck
R CMD check --timings --as-cran MARSS_3.11.7.tar.gz
```

## Check on other OS

* [ ] `devtools::check_win_devel()`
* [ ] `rhub::check_for_cran()`


## Add cran-comments

* [ ] Update `cran-comments.md`


## Prepare for release:

* [ ] Check [current CRAN check results](https://cran.rstudio.org/web/checks/check_results_MARSS.html)
* [ ] [Polish NEWS](https://style.tidyverse.org/news.html#news-release)
* [ ] [`urlchecker::url_check()`](https://github.com/r-lib/urlchecker)


Re-read the submission guidelines to see if anything changed: https://cran.r-project.org/submit.html

Check that package builds on Windows

- [ ] Upload the Dropbox/MARSS_X.X.X.tar.gz file to https://win-builder.r-project.org

Upload package to CRAN

- [ ] https://cran.r-project.org/submit.html

