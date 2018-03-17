MARSS <img src="logot.png" align="right" />
=====================================================

[![cran version](http://www.r-pkg.org/badges/version/MARSS)](https://cran.r-project.org/package=MARSS)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/MARSS?)](https://github.com/metacran/cranlogs.app)
[![](https://cranlogs.r-pkg.org/badges/MARSS?color=FFD700)](https://www.r-pkg.org/pkg/MARSS)

MARSS stands for Multivariate Auto-Regressive(1) State-Space. The MARSS package is an R package for estimating the parameters of linear MARSS models with Gaussian errors.  This class of model is extremely important in the study of linear stochastic dynamical systems, and these models are important in many different fields, including economics, engineering, genetics, physics and ecology.  The model class has different names in different fields, for example in some fields they are termed dynamic linear models (DLMs) or vector autoregressive (VAR) state-space models.  The MARSS package allows you to easily fit time-varying constrained and unconstrained MARSS models with or without covariates to multivariate time-series data via maximum-likelihood using primarily an EM algorithm.

### Collaborate

Issues? [https://github.com/nwfsc-timeseries/MARSS/issues]()
Wiki [https://github.com/nwfsc-timeseries/MARSS/wiki]()

### INSTALL

To install MARSS from CRAN:

```
install.packages("MARSS")
library(MARSS)
```

The latest release on GitHub may be ahead of the CRAN release. To install the latest release on GitHub:
```
install.packages("devtools")
library(devtools)
install_github("nwfsc-timeseries/MARSS@*release")
library(MARSS)
```

The master branch on GitHub has work leading up to a GitHub release.  The code here may be broken though usually prelim work is done on a development branch before merging.  To install the master branch:
```
install_github("nwfsc-timeseries/MARSS")
```


### DOCUMENTATION and TUTORIALS

There is an extensive user manual included in the package:
[here](https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf).

Many applications are also covered in our Applied Time Series Analysis book: [here](https://nwfsc-timeseries.github.io/atsa-labs/).

We have lectures on our course website: [here](https://nwfsc-timeseries.github.io/atsa/).


### CITING:

If you use MARSS results in publications, please cite the primary citation:

Holmes, E. E., Ward, E. J. and Wills, K. (2012) MARSS: Multivariate Autoregressive State-space Models for Analyzing Time-series Data. The R Journal. 4(1):11-19

You can also cite the package as you would other R packages:

Elizabeth Holmes, Eric Ward, Mark Scheuerell, and Kellie Wills (2018). MARSS: Multivariate Autoregressive State-Space Modeling. R package version 3.10.4.

Update the version number and year if you use a more recent version on GitHub.

### PUBLICATIONS

To see our publications using MARSS models, see the [NWFSC Time-Series Analysis website](https://nwfsc-timeseries.github.io/).