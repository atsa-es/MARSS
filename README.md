MARSS
=============
[![cran version](http://www.r-pkg.org/badges/version/MARSS)](https://cran.r-project.org/package=MARSS)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/MARSS?)](https://github.com/metacran/cranlogs.app)

MARSS stands for Multivariate Auto-Regressive(1) State-Space. The MARSS package is an R package for estimating the parameters of linear MARSS models with Gaussian errors.  This class of model is extremely important in the study of linear stochastic dynamical systems, and these models are important in many different fields, including economics, engineering, genetics, physics and ecology.  The model class has different names in different fields, for example in some fields they are termed dynamic linear models (DLMs) or vector autoregressive (VAR) state-space models.  The MARSS package allows you to easily fit time-varying constrained and unconstrained MARSS models with or without covariates to multivariate time-series data via maximum-likelihood using primarily an EM algorithm.

### INSTALL

To install MARSS from CRAN:

```
install.packages("MARSS")
library(MARSS)
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

Elizabeth Holmes, Eric Ward, and Kellie Wills (2013). MARSS: Multivariate Autoregressive State-Space Modeling. R package version 3.9.

Update the version number and year if you use the more recent version on GitHub.

### PUBLICATIONS

To see our publications using MARSS models, see the [NWFSC Time-Series Analysis website](https://nwfsc-timeseries.github.io/).