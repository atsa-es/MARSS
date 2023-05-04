<style>
.nav{
    border:1px solid #ccc;
    border-width:1px 0;
    list-style:none;
    margin:0;
    padding:0;
    text-align:center;
}
.nav ul{
    margin-before: 1em;
    margin-after: 1em;
    margin-start: 0;
    margin-end: 0;
    padding:0;
}
.nav li{
    display:inline-block;
}
.nav a{
    display:inline-block;
    padding:5px;
}
</style>

<ul class="nav">
  <li><a href="#install">Install</a></li>
  <li><a href="#documentation">Documentation</a></li>
  <li><a href="https://github.com/atsa-es/MARSS/issues">Issues</a></li>
  <li><a href="#cite">Citation</a></li>
  <li><a href="#pubs">Publications</a></li>
  <li><a href="#license">License</a></li>
  <li><a href="https://atsa-es.github.io/MARSS/NEWS.html">NEWS</a></li>
</ul>

<img src='man/figures/logo.png' align="right" height="139" style="margin:15px 10px"/>

MARSS stands for Multivariate Auto-Regressive(1) State-Space. The MARSS package is an R package for estimating the parameters of linear MARSS models with Gaussian errors.  This class of model is extremely important in the study of linear stochastic dynamical systems, and these models are important in many different fields, including economics, engineering, genetics, physics and ecology.  The model class has different names in different fields, for example in some fields they are termed dynamic linear models (DLMs) or vector autoregressive (VAR) state-space models.  The MARSS package allows you to easily fit time-varying constrained and unconstrained MARSS models with or without covariates to multivariate time-series data via maximum-likelihood using primarily an EM algorithm.

[![cran version](http://www.r-pkg.org/badges/version/MARSS)](https://cran.r-project.org/package=MARSS)
[![github](https://img.shields.io/github/v/release/atsa-es/MARSS?color=brightgreen&label=GitHub)](https://github.com/atsa-es/MARSS/releases/latest)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/MARSS?)](https://github.com/metacran/cranlogs.app)
[![](https://cranlogs.r-pkg.org/badges/MARSS?color=FFD700)](https://www.r-pkg.org/pkg/MARSS)
[![License: CC0-1.0](https://img.shields.io/badge/License-CC0%201.0-lightgrey.svg)](http://creativecommons.org/publicdomain/zero/1.0/)


### INSTALL {#install}

To install MARSS from CRAN:

```
install.packages("MARSS")
library(MARSS)
```

The latest release on GitHub may be ahead of the CRAN release. To install the latest release on GitHub:
```
install.packages("devtools") # if needed
devtools::install_github("atsa-es/MARSS@*release")
library(MARSS)
```

The master branch on GitHub has work leading up to a GitHub release.  The code here may be broken though usually prelim work is done on a development branch before merging.  To install the master branch:
```
install_github("atsa-es/MARSS")
```
If you are on a Windows machine and get an error saying 'loading failed for i386' or similar, then try
```
options(devtools.install.args = "--no-multiarch")
```
If R asks you to update packages, and then proceeds to fail at installation because of a warning that a package was built under a later R version than you have on your computer, use
```
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE)
```
To install an R package from GitHub, you need to be able to build an R package on your machine. If you are on Windows, that means you may need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/). In more recent versions of R, it seems like the Rtools dependency for Windows users has been removed, so try installing. If you get an error about no gcc installation, it means you need Rtools. On a Mac, installation should work fine; you do not need to install anything.

### DOCUMENTATION and TUTORIALS  {#documentation}

- [Quick Start Guide](https://CRAN.R-project.org/package=MARSS/vignettes/Quick_Start.pdf).
- [User Manual](https://CRAN.R-project.org/package=MARSS/vignettes/UserGuide.pdf) - The extensive user manual included in the package has many examples of how to fit MARSS models to a variety of data sets.
- [ATSA lab book](https://atsa-es.github.io/atsa-labs/) - Many applications are also covered in our Applied Time Series Analysis book developed from the labs in our course.
- [ATSA course website](https://atsa-es.github.io/atsa/) - We have lectures and all material from our course on our course website. Select the Lectures tab to find the lecture material and videos of lectures.
- [Wiki](https://atsa-es.github.io/MARSS-wiki/) - The MARSS wiki has misc and random projects and code.


### ISSUES and BUG REPORTS {#bugs}

Issues? [https://github.com/atsa-es/MARSS/issues](https://github.com/atsa-es/MARSS/issues)

### CITATION  {#cite}

If you use MARSS results in publications, please cite the primary citation:

Holmes, E. E., Ward, E. J. and Wills, K. (2012) MARSS: Multivariate Autoregressive State-space Models for Analyzing Time-series Data. The R Journal. 4(1):11-19

You can also cite the package and user guide:

Elizabeth E. Holmes, Eric J. Ward, Mark D. Scheuerell and Kellie Wills (2020). MARSS:
  Multivariate Autoregressive State-Space Modeling. R package version 3.11.4.
  
Holmes, E. E., M. D. Scheuerell, and E. J. Ward (", year, ") Analysis of multivariate time-series using the MARSS package. Version ", meta$Version, ". NOAA Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112, DOI: 10.5281/zenodo.5781847
  
Type `citation("MARSS")` at the command line to get the most up to data citations.

### PUBLICATIONS {#pubs}

To see our publications using MARSS models, see the [Applied Time Series Analysis website](https://atsa-es.github.io/).


### NOAA Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and
Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is
provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of
Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or
favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a
DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by
DOC or the United States Government.
