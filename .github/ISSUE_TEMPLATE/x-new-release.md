---
name: New Release
about: Steps for a new MARSS release to CRAN
title: ''
labels: ''
assignees: eeholmes

---

- [ ] Run checks without building the vignettes
```
devtools::check(document = FALSE, vignettes=FALSE)
```
- [ ] Run examples only
```
devtools::run_examples(pkg="~/Documents/GitHub/MARSS", document=FALSE, fresh=TRUE)
```

- [ ] spellcheck the whole package
```
#devtools::install_github("ropensci/spelling")
spelling::spell_check_package()
fils <- dir("vignettes/manual_files", full.names = TRUE)
fils <- fils[stringr::str_detect(fils, "Rnw")]
spelling::spell_check_files(fils, ignore = spelling::get_wordlist(pkg = "."), lang = "en_US")
```

- [ ] Check Manual
```
setwd("~/Documents/GitHub")
devtools::build_manual(pkg="MARSS")
```

- [ ] Build the package into tar.gz. Note do not use `devtools::build()` it destroys the `inst/doc` folder.
Note, the Mac LaTeX installation doesn't install packages on the fly. You may need to install the needed packages. From terminal: `tlmgr install imakeindex`, `tlmgr install collection-fontsrecommended`, `tlmgr install footmisc`, tlmgr install appendix`, plus a few more in header.tex.
 
Open terminal
```
cd ~/Documents/GitHub
rm MARSS_3.11.4.tar.gz
R CMD build MARSS
rm -r ~/Dropbox/MARSS
rm ~/Dropbox/MARSS_3.11.4.tar.gz
cp MARSS_3.11.4.tar.gz ~/Dropbox/MARSS_3.11.4.tar.gz
cd ~/Dropbox
tar -xvzf MARSS_3.11.4.tar.gz
rm ~/Dropbox/MARSS/vignettes/Makefile
R CMD build --no-build-vignettes MARSS
```

Run version tests
- [ ] Run this `vignettes/versiontest.R`

- [ ] Run checks on the built package
```
rm -r ~/Dropbox/MARSS.Rcheck
R CMD check --timings --as-cran MARSS_3.11.4.tar.gz
```

Run internal tests 
  - [ ] Run `tests/model.R` to set up models for tests
  - [ ] Then run `devtools::test()` but I find running the tests one by one in RStudio is easier to make sense of. Some tests will through errors, but all test should pass. Errors are correct if the test passes and result should be try-error.
       - [ ] test.coef
       - [ ] test-forecast
       - [ ] test-GDF
       - [ ] test-KFAS
       - [ ] test-MARSSkfss
       - [ ] test-plotting
       - [ ] test-print
       - [ ] test-residualsHarvey
       - [ ] test-structTS
       - [ ] test-tt

Check the documentation pdfs
  - [ ] EMDerivation
  - [ ] Quick Start
  - [ ] Residuals
  - [ ] User Guide - See separate issue for this

Compress documentation pdfs
  - [ ] Manually compress the pdfs in the doc folder Dropbox/MARSS/inst/doc
  - [ ] Rebuild package with smaller docs, and then re-Run checks on the built package to see that it passes size checks
```
R CMD build --no-build-vignettes MARSS
rm -r ~/Dropbox/MARSS.Rcheck
R CMD check --timings --as-cran MARSS_3.11.4.tar.gz
```

Update the doc folder for GitHub. This has to be done manually by copying from
  - [ ] Dropbox/MARSS/inst/doc into GitHub/MARSS/inst/doc

Transition MARSS to atsa-es
  - [ ] Clone v3.11.4 over to atsa-es so that url links pass. Leave MARSS on nwfsc-timeseries for the time being.

Prepare for release:
* [ ] Check [current CRAN check results](https://cran.rstudio.org/web/checks/check_results_MARSS.html)
* [ ] [Polish NEWS](https://style.tidyverse.org/news.html#news-release)
* [ ] [`urlchecker::url_check()`](https://github.com/r-lib/urlchecker)

Update the doc folder on GitHub. This has to be done manually
- [ ] Copy Dropbox/MARSS/inst/doc into GitHub/MARSS/inst/doc

Re-read the submission guidelines to see if anything changed: https://cran.r-project.org/submit.html

Check that package builds on Windows
- [ ] Upload the Dropbox/MARSS_X.X.X.tar.gz file to https://win-builder.r-project.org

Upload package to CRAN
- [ ] https://cran.r-project.org/submit.html



The items below are for cases where you are submitting the package files without any extra processing. I process my vignette PDFs separately, so the below won't work. But when submitting there is a link to check the tar.gz file on Windows.



* `devtools::check_win_devel()`
* `rhub::check_for_cran()`
*  `revdepcheck::revdep_check(num_workers = 4)`

Submit to CRAN:

* `usethis::use_version('patch')`
* `devtools::submit_cran()`
* Approve email

Wait for CRAN...

* [ ] Accepted :tada:
