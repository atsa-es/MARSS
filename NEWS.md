---
output: html_document
---

MARSS Development site
------------------------------------
New work on MARSS before posting to CRAN is at the GitHub repo.  See issues posted there.

Status
-----------------------------------

7-09-2020

* The versiontest.R tests passed.
* Working on testing predict and residuals against other packages and models.
* Currently working on `StructTS()` examples in Chapter_Structural_TS.Rnw. 

7-12-2020

* Working on issue #74 that is coming up in the BSM model. Weird Q update error at iteration 1. It's a numerical accuracy issue. Fixed by using solve(A,B) to get J0 instead of inverse of Vtt1[,,1] [Fixed]
* Rerun the versiontest.R tests because of the J0 change. [passed 7/14]

7-14-2020

* Back to working on testing predict and residuals against other packages and models. After fixing issue #74 and ensuring that version test passes.
* Currently working on `StructTS()` examples in Chapter_Structural_TS.Rnw. 

7-16-2020

* Fixed various issues in MARSSresiduals() due to not passing type in. That error progated to problems in plot functions. 
* autoplot.marssPredict() not working for forecasts since facet_wrap fails when plot uses data with different number of time steps. Fixed by not subsetting but instead using NAs for the data I don't want to show.
* Did the forecast subsection for StructTS models.

7-20-2020

* Fixed bug in MARSSresiduals.tT. Needed t(chol()) for Cholesky standardization.
* Working on innovations state residuals. Finished. Residuals.Rnw and MARSSresiduals.tt1()
* Changed ACF plots to be for the innovations residuals. Smoothation residuals are temporally correlated.

7-22-2020

* Revamped fitted, tidy and predict. 
* fitted.marssMLE -> MARSSpredict. 
* tidy.marssMLE -> fitted.marssMLE (the observations and states bit). tidy only return parameter info
* predict.marssMLE updated to use ytt1, ytt, ytt1. Only minor changes to predict.marssMLE.

To do

* Work on fitted() and residuals() for StructTS models.
* Then move to multivariate examples.
* KFAS examples: https://www.rdocumentation.org/packages/KFAS/versions/1.3.7/topics/KFAS

7-23-2020

* Back to original. Got rid of MARSSpredict and MARSSest
* use tsStruct for the smoothed and filter ests

* need to fix equations in Rd files. changed to x(t+1) and messed up time indexing for the 
time-varying parameters.
x(t) = B(t)*x(t-1)+u(t)+C(t)c(t)
so 
x(t+1) = B(t+1)x(t)+u(t+1)
Go back to x(t). That's how it is in the User Guide.

* In Covariates.Rnw I show acf of residuals. need to use innovations there.
* In quick examples, I talk about diagnostics too.
* changed default residuals(fit) to return innovations.
    * don't return the state innovation residuals (only via MARSSresiduals())
    * check all refs to residuals in the documentation [done but recheck text]
* Check getDFAfits(). does it return smoothations like the chap says it does?

    
MARSS 3.11.0 (resids_update for CRAN)
------------------------------------
Version 3.11.0 is focused on the `predict`, `fitted` and `residuals` functions functions. Most of the `predict` changes are listed below for 3.10.13 release on GitHub.


ENHANCEMENTS

* `ldiag()` convenience function added to make list diagonal matrices. This replaces having to do code like `a <- matrix(list(0),2,2); diag(a) <- list(2,"a")`. Now you can call `ldiag(list(2,"a"))`.
* Added `accurancy.marssMLE()` and `accuracy.marssPredict()` which returns accuracy metrics sensu the **forecast** package.
* Added `is.unitcircle()` utility function and added tol so that it does not fail if abs(eigenvalue) is above 1 by machine tolerance.
* Added ACF plots for model and state residuals to `plot.marssMLE()` and `autoplot.marssMLE()`.
* Revamped `residuals.marssMLE()`. Got rid of `augment.marssMLE()` and renamed it `residuals.marssMLE()`. The old `residuals.marssMLE()` became `MARSSresiduals()`. There was too much duplication between `residuals.marssMLE()` and `augment.marssMLE()` and between `augment.marssMLE()` and `fitted.marssMLE()`. Also I want to minimize dependency on other packages and `augment` is a class in the **broom** package. This required changes to the `glance.marssMLE()`, `plot.marssMLE()` and `autoplot.marssMLE()` code.
* Revamped `fitted.marssMLE` and `tidy.marssMLE`. `fitted.marssMLE` now returned the estimates from the Kalman filter or smoother which `tidy.marssMLE` had returned. `tidy.marssMLE` only returns a data frame for the parameter estimates.
* V0T was computed with an inverse of Vtt1[,,1]. This led to unstable numerics when V00 was like matrix(big, m, m). Changed to use `solve(t(Vtt1[,,1]), B%*%V00)` which should be faster and seems to have lower numerical error.

BUGS

* This bug affected `residuals()` in cases where R=0. In v 3.10.12, I introduced a bug into `MARSSkfss()` for cases where R has 0s on diagonal. **History**: To limit propogation of numerical errors when R=0, the row/col of Vtt for the fully determined x need to be set to 0. In v 3.10.11 and earlier, my algorithm for finding these x was not robust and zero-d out Vtt row/cols when it should not have if Z was under-determined. This bug (in < 3.10.12) only affected underdetermined models (such as models with a stochastic trend and AR-1 errors). To fix I added a utility function `fully.spec.x()`. This returns the x that are fully determined by the data. There was a bug in these corrections which made `MARSSkfss()$xtT` wrong whenever there were 0s on diagonal of R. This would show up in `residuals()` since that was using `MARSSkfss()` (in order to get some output that `MARSSkfas()` doesn't provide.) The problem was `fully.spec.x()`. It did not recognize when Z.R0 (the Z for the R=0) was all 0 for an x and thus was not (could not be) fully specified by the data. Fix was simple check that colSums of Z.R0 was not all 0.
* When computing the Cholesky standardized residuals, the lower triangle of the Cholesky decomposition should be used so that the residuals are standardized to a variance of 1. `base::chol()` returns the upper triangle. Thus the lines in `MARSSresiduals.tT()` and `MARSSresiduals.tt1()` that applied the standardization need `t(chol())`. 
* `trace=1` would fail because `MARSSapplynames()` did not recognize that `kf$xtt` and `kf$Innov` were msg instead of a matrix. I had changed the `MARSSkfas()` behavior to not return these due to some questions about the values returned by the KFAS function.
* `is.validvarcov()` used eigenvalues >= 0 as passing positive-definite test. Should be strictly positive so > 0.
* `MARSSkfas()` had bug on the line where `V0T` was computed when `tinitx=0`. It was using `*` instead of `%*%` for the last `J0` multiplication. It would affect models with a non-zero `V0` under certain `B` matrices, such as structural models fit by `StructTS()`.
* `MARSSresiduals.tT()` and `MARSSresiduals.tt1()` but was in old `residuals.marssMLE()` also. If MLE object had the `kf` element, `kf` was not assigned in code since there was no `kf <- MLEobj$kf` line for that case. Normally MLE objects do not have the `kf` element, but it could be added or is added for some settings of `control$trace`.

DOCUMENTATION and MAN FILES

* Added covariates and example to `MARSS_dfa.Rd`
* Added chapter on Structural time series models which compares `StructTS()` to `MARSS()` output.
* Removed all mention of `augment()` from documentation and manuals. Replaced with `residuals()`.
* `predict.marssMLE.Rd` (help page) had bug in the examples. remove `Q=Q` from the model list in the first example.
* Cleaned-up the man pages for `predict()` and `predict.marssMLE()`.
* In the chapter on structural breaks and outliers, Koopman et al (1998) use the marginal residuals in their example rather than the Cholesky standardized residuals. Changed to use marginal residuals to follow their example.

OTHER

* `fitted.marssMLE` now called `MARSSpredict`. It is the expected value of the right side of the MARSS equations.
* The part of `tidy.marssMLE` related to states and observations was moved to `fitted.marssMLE`. It is the expected value of the left side of the MARSS equations.
* Changed `fitted.marssMLE` to have column with xtt when type="xtt1" instead of xtT.
* Changed the x0 estimation behavior for `predict.marssMLE()` when no data passed in.
* Added x0 argument to `predict.marssMLE()` so that user can specify x0 if needed.
* Removed the tibble class from the data frames returned by `residuals.marssMLE()`. The data frames are still in tibble form. Removed all reference to tibbles in the documentation.
* When `trace = -1` some tests were still being done. I added a test for `trace = -1` to a few more test lines in `MARSS.R` and `MARSS_marxss.R`.


MARSS 3.10.13 (GitHub 2-25-2020)
------------------------------------
Version 3.10.13 mainly has to do with the `predict()` and `forecast()` functions along with plotting and printing methods.

ENHANCEMENTS

* `MARSSkfss()` and `MARSSkfas()` Add rownames to the x elements of the list.
* `MARSSkf()` Added `newdata` to allow user to pass in a new dataset to fit with the fitted model.
* `predict.marssMLE()` Shows the prediction or confidence intervals for data or states. Forecasts can be done by passing in `h`. `newdata` can be passed in also and fitted model will be used to fit these data and show the intervals. Output is same form as a tibble (but not a tibble). Returns a list of class `marssPredict`.
* `forecast.marssMLE()` This does the foreward forecasting past the end of the data. Intended to be called by `predict.marssMLE`. I did not write a `marssMLE` method for the `forecast` generic in the **forecast** package since that would require that the forecast package be required for the **MARSS** package.
* `plot.marssPredict()` plot method for the new marssPredict object. This is designed to look like the `plot.forecast()` function in the **forecast** package.
* `print.marssPredict()` print method for marssPredict objects.

MARSS 3.10.12 (CRAN 2-3-2020)
------------------------------------
Version 3.10.12 update mainly has to do with the `tidy()`, `fitted()` and `augment()` enhancements which clarify ytT, xtT and residual intervals for MARSS models. This was a major update though probably users will not notice much as it only affects residuals output. A few minor bugs were fixed which caused errors to be thrown in some rare time-varying cases. One bug that affected bootstrap confidence intervals was fixed. The documentation got a major clean-up. The Residuals report has been heavily edited to improve precision and clarity (with added verbosity). The help files and automated manual from the help files were cleaned-up  and  some of the internal functions moved out of the manual.

BUGS

* `MARSSsimulate()` missing values were placed in the wrong positions in simulated data. This would affect all simulated data with missing values and thus any function that used `MARSSboot()`, for example bootstrap confidence intervals for a data set with missing values. Default is to use Hessian so the user would not normally have encountered this bug and it had little effect on the CIs.
* `fitted.marssMLE()` Fixed bug in fitted.marssMLE for states when one.step.ahead=TRUE. It was using xtt1[,t-1] instead of xtt[,t-1]. The former meant it only used data up to t-2.
* `degen.test()` in MARSSkem() was not catching when R or Q was time-varying (and thus degeneracy should not be allowed). Changed to test the 3D of model.dims == 1 or not.
* `residuals.marssMLE(..., Harvey=TRUE)` would fail if Q, B, or G was time-varying because parmat() called with t+1. Changed to only call parmat() when t<TT. Q, B, and G do not appear in the recursion when t=TT so parmat() with t=t+1 is never needed.
* `MARSS_dfa()` used form="dfa" in MARSS.call list. Just info. Never used.
* Default A matrix ("scaling") was throwing an error for manually set up DLM models. Problem was call to check that Z was a design matrix in MARSS_marxss.R. It was not catching that Z was time-varying before running `is.design()`.
* `toLatex.marssMODEL()` Fixed some old bugs in toLatex_marssMODEL.R. Added S3 class declaration in NAMESPACE for toLatex. fixed equation attribute in MARSS_marxss. G{t} was used instead of G_{t}. Only affected toLatex_marssMODEL(). Had extra line in build.mat.tex() that removed last line of matrices. This function was not exported so users would never have run into these bugs.
* `MARSSkfss()` To limit propogation of numerical errors when R=0, the row/col of Vtt for the fully determined x (determined from data) need to be set to 0. My algorithm for finding these x was not robust and zero-d out Vtt row/cols when it should not have if Z was under-determined. MARSSkfss() is not used for fitting and only affected underdetermined models (such as models with a stochastic trend and AR-1 errors). To fix I added a function `fully.det.x()` to the utility functions. This returns the x that are fully determined by the data. Note, MARSkfss() is the classic Kalman filter/smoother. The MARSS algorithm does not use this normally. Normally MARSSkfas(), build off the Koopman et al algorithm which avoids unneeded matrix inverses, is used. MARSSkfas() uses the Kalman filter/smoother in the KFAS package. 

ENHANCEMENTS

* `MARSShatyt()` Added ytt, ytt1, Ott, Ott1 to MARSShatyt() so that tidy.marssMLE() can more easily return the one-step-ahead preditions for Y(t). Also added var.ytT and var.EytT so you can easily get the estimates, CI and prediction intervals for missing data. Added only.kem to MARSShatyt() so that only values conditioned on 1:T as needed by MARSS kem are returned. This makes the Ey part of a MARSS object smaller and speeds up MARSShatyt() a little.
* `tidy.marssMLE()` Changed type for tidy() to xtT, ytT and fitted.ytT. tidy() exclusively gives estimates of things (parameters, X, Y, fitted Y) conditioned on all the data.
* `fitted.marssMLE()` Added interval=c("none", "confidence", "prediction") to fitted() and returns a list with se's (or sd's if prediction) and intervals. Also added conditioning argument to fitted.marssMLE which gives fitted values with different conditioning. Changed default output to tibble.
* `augment.marssMLE()` Changed standard errors output for augment() to .se.fit for std error of fitted y and .sigma to std error of residuals. This matches what augment.lm outputs.
* `plot.marssMLE()` and `autoplot.marssMLE()`. Added plot.par to plot.marssMLE and autoplot.marssMLE so that the plots can be customized. Added plot of ytT to both functions. Changed the residuals plots to use the CIs for the residuals not the loess CIs.
*  `MARSSinfo()` Added "AZR0" to MARSSinfo() to give info if user gets error that A cannot be estimated with R=0.  Added more informative message to MARSSkemcheck() for that case.
* Made all if statements checking class of object robust to the class returning more than one class (so vector of length > 1). Due to change in R 4.0.0 where matrix has class c("matrix","array")
* Updated all code to tidyverse style
* Changed `residuals.marssMLE()`. This is now a helper function which calls `MARSSresiduals.tT()` or `MARSSresiduals.tt1()`. The former is smoothation residuals and the latter is innovations (one-step-ahead) residuals.

DOCUMENTATION and MAN FILES

* Added derivation of variance of Y conditioned on y and X to EMDerivation.Rnw.  Needed for CI on the missing values estimate.
* Major update to Residuals report. No changes to equations but much editting to improve precision and clarity (with much more verbosity). Reposted to Arxiv. Added innovations residuals.
* tidy, augment and fitted man files got major update.
* internal functions given \keyword{internal} so they don't appear in the documentation, but will appear if you use `?` or help.search.
* Rd files extensively cleaned to improve linking and move more of the internal functions out of view of the normal user.  Equations cleaned up (though not completely).

MARSS 3.10.11 (GitHub 8-3-2019)
------------------------------------
Minor update. Version 3.10.11 has some edits to speed up the code by minimizing calls to expensive checking functions and fixes a bug in `MARSSharveyobsFI()` that appeared if a parameter was fixed and time-varying and `MARSSparamCIs()` was called.

BUGS

* `MARSSkfss()` had a bug "\*" was used in place of %\*% with J0. Would never show up unless V0 estimated.
* Bug in `MARSSharveyobsFI()` which arose if a parameter was fixed and time-varying. This caused MARSSparamCIs() to fail if a parameter was fixed and time-varying.
* `dparmat()` did not return values if it was time-varying and fixed. Caused tidy() to return error for dlm models.

ENHANCEMENTS

* `is.validvarcov()` is expensive. minimize calls to it. If called with a diagonal matrix, it should automatically pass so added a check to is.validvarcov() to see if matrix is diagonal.
* `is.marssMLE()` is expensive. Replace with a call to class().
* Added S3 methods for broom functions if broom is loaded.
* Added `autoplot.marssMLE()` function and updated plot documentation to cover autoplot functions.
* Fixed typos in Case Study 4 and Derivation eq 143b.


MARSS 3.10.10 (CRAN 11-2-2018)
------------------------------------
Minor update to declare S3 objects if user has broom package installed.  A few minor changes also made.

ENHANCEMENTS

* Added if statement to declare S3 objects in NAMESPACE if user has broom package installed.  Allows delayed loading of S3 methods for broom.
* `is.validvarcov()` is expensive.  Minimize use in MARSSaic, MARSSboot, MARSSoptim, MARSSparamCIs, MARSSsimulate, MARSS_marxss.
* Added `is.diagonal()` utility function
* Typo in User Guide Chapter 2.  R was supposed to be "diagonal and equal" in the first example.
* Changed `plot.marssMLE()` to `autoplot.marssMLE()` since it is ggplot2 based
* Created `plot.marssMLE()` that uses only base R graphics

BUGS

* `MARSSkfss()` had bug on V0T line.  "\*" instead of %\*%.  This code is almost never called.
* `MARSSkfss()` had bug.  Did not check if G was time-varying.  This code is almost never called.


MARSS 3.10.8 (CRAN 4-14-2018)
------------------------------------
Major update over 3.9. The main changes have to do with with errors in the Hessian matrix whenever the Cholesky of the R or Q matrix was used (when they weren't diagonal).  This affected all the residuals and confidence intervals calculations for non-diagonal R and Q. Hessian for non-diagonal Z was also bad. Version 3.10.8 completely abandons working with the chol transformed variance-covariance matrices for the Hessian calculation.  The chol transformation was not necessary for computing the Hessian since the Hessian is computed at the MLEs and localized. Also the default Hessian computation now uses the Harvey et al. analytical algorithm for the Hessian rather than a numerical estimate.

BUGS

`residuals.marssMLE()`

* Erroneous standardized residuals when Z was non-diagonal (and thus also > 1 row). 
* Changed residuals.MARSSMLE so it returns residuals (which would be equal to 0) when there are 0s on the diagonal of Q or R
* If t=1 and x0 was fixed, residuals were not returned.
* The wrong residuals were returned if there were any missing values.  Fix involved implementing the missing values modifications for the Kalman filter described in Shumway and Stoffer.

inits functions

* In `MARSSinits_marxss()` function would give error if U, A, C, or D fixed and user passed in inits.  inits ignored in this case so should not throw error.
* `alldefaults` could be updated by form.  A few functions were neglecting to (re)load alldefaults or to ressign `alldefaults` when updated: is_marssMLE(), MARSSinits.marxss(), MARSSinits().  The variables in the pkg_globals environment should be (and only be) loaded when needed by a function and only loaded into the function environment.

kalman filter functions

* `MARSSkf()` was not passing optional function args to `MARSSkfas()`.
* `MARSSkfss()` mis-counting num data points when R=0, V0=0, and tinitx=1.  When Ft[,,1]=0 (e.g. when R=0, V0=0, and tinitx=1), MARSSkfss() was including the y[1] associated with Ft[,,1]=0 in the # number of data points.  These should be excluded since they don't affect x10.

Confidence intervals and std error for R and Q

* `MARSSparamCIs()` gave the wrong s.e. for variances and covariances when method="hessian".  It also gave the wrong CIs for variances and covariances when var-cov matrix was non-diagonal. There were a series of issues related to back-transforming from the hessian of a chol-transformed var-cov matrix ( Sigma=chol%\*%t(chol) ).
* `MARSSparamCIs()`, vrs 3.9 was getting the Hessian matrix numerically using a var-cov matrix that had been transformed with a Cholesky decomposition to ensure it stays positive-definite.  The upper and lower CIs were computed from the s.e.'s. I back-transformed the Hessian to the original (non-chol) scale the same way I back transformed a var-covariance matrix.  But the variance of s^2 is not the var(s)^2, which is what I was doing, essentially. So the s.e. for R and Q were wrong in all cases.  Note, using a Hessian to estimate CIs for variance-covariance matrices is generally a bad idea anyhow however.
* For non-diagonal matrices.  There was a bug in `MARSShessian()` in subscripting the d matrix when doing the chol transform.  Caused NAs for cases with non-diagonal matrices.  However, had the se been returned, they would be wrong for non-diagonal matrices because the elements of the chol transformed matrices do not correspond one-to-one to the non-transformed matrices.  E.g. the untransformed [2,2]^2 is chol transformed [1,2]^2+[2,2]^2.  The hessian used was for the chol-transformation and the curvature of the LL surface for the chol-transformed values is different than the curvature for the untransformed var-cov elements.

Fix: I completely abandoned working with the chol transformed variance-covariance matrices for the Hessian calculation.  The chol transformation was not necessary for computing the Hessian since the Hessian is computed at the MLEs and localized.  

1. Created a new function `MARSSharveyobsFI()` which uses the Harvey (1989) recursion to analytically compute the observed Fisher Information matrix. This is the Hessian for the untransformed var-cov parameters.  So CIs on variances can be negative since the variance of the MLE is being approximated by a MVN (which can lead to negative lower CIs).

2. Harvey1989 is now the default function when method='hessian'. *In later vrs of MARSS, this is changed to Holmes2014.*

3. The user can also select method='hessian' and hessian.fun='fdHess' or hessian.fun='optim'.  This will compute the Hessian (of the log-LL function at the MLEs) numerically using these functions.  The variance-covariance matrices are NOT chol transformed.  These are numerically estimated Hessians of the untransformed variance-covariance matrices.

4. Added `MARSSinfo(26)` which discusses the reason for NAs in the Hessian.

MISC MINOR BUGS

* `MARSShatyt()` was setting ytt1 (expected value of y(t) conditioned on data up to t-1) to y(t), which is incorrect. Expected value of y(t) conditioned on y up to t-1 is Z xtt1 + a
* `CSEGriskfigure()` panel 2 was wrong when mu>0 (increasing population).
* `CSEGriskfigure()` panel 2 CIs were wrong when CI method=hessian since Q had not been back transformed (so was using sqrt(Q)).
* `MARSSkemcheck()`'s test that fixed B is in unit circle failed when B was time-varying and some fixed and others estimated. Also when the eigenvalues were complex, it should test the real part only.
* `coef.marssMLE()` was not stopping when illegal "what: arg passed in.  man page did not say what happens when type=parameter.
* Passing in method not in allowed.methods was causing errors since model conversion and testing happening before checkMARSSinputs and model testing is algorithm dependent (some forms not allowed by BFGS).  Added check for method at top of MARSS function.

IMPROVEMENTS

* w(t) and v(t) can be specified as G(t)\*w(t) and H(t)\*v(t) where G and H are fixed matrices (not estimated).  In version 3.10, G and H are restricted to being 0 or identity, however the code is in place for other values.
  - changes to `MARSS.marxss()` and `MARSS.marss()` to allow G, H, and L passed in
  - change to `MARSSkem()` to specifiy star lists with G, H, and L (mathbb(elem) in EM Derivation)
  - changed `MARSSkss()` to use Q\*=G Q t(G), R\*=H R t(H) and V0\*=L V0 t(L)
* Removed the function `MARSSmcinits()` and added chapter on searching over the initial conditions into the User Guide.  As the MARSS models that MARSS() can fit expanded, `MARSSmcinits()` was increasing obsolete and it was impossible to come up with good searching distributions.  Because `MARSSmcinits()` was removed, control$MCInits list item was removed also from defaults and from accepted input.
* Added default inits for c and d in marxss form so that user can pass in inits using coef(fit); was balking because this includes d and c which didn't have defaults.  Removed msg referring to need that model be in marss form for inits (not true).
* Changed `MARSS.marxss()` to allow c and d to be 3D arrays.  This allows one to use inits=fit to set inits and not get a d (or c) must be 2D error.
* Added info to `MARSSinfo(4)` regarding errors about R=0 and x0 not fixed. Added info to the error warnings to direct user to MARSSinfo().
* Changed order of MARSS args to be MARSS(y, model= , inits=, ...)
* Added `pchol()` and `psolve()` functions to return the chol or inverse via solve when there are 0s on the diagonal
* Added information to print.marssMODEL on summary.marssMODEL.  Added silent argument to summary.marssMODEL to block printing to the console.
* `print.marssMLE(x, what="par")` returned a vector of estimated values instead of the list of par.  Changed to return the list.
* Added E[y(t), x(t+1)] to `MARSShatyt()` output.  Needed for residuals.marssMLE().
* Added code to `is.validvarcov()` so that it returns an error if the user specifies a structurally illegal variance-covariance matrix.  Added info to MARSSinfo(25).
* Added a model.frame method for marssMODEL and marssMLE
* Added broom `augment`, `tidy` and `glance` functions for marssMLE
* Added `logLik` method for marssMLE objects
* Added `fitted` method for marssMLE objects to return Z xtT + u (model fitted value of y)
* Added `plot` method with diagnostics

DOCUMENTATION

* Added Multivariate linear regression chapter
* Added chapter on estimating a Leslie matrix from stage time series using a MARSS model.
* Added chapter on searching over the initial conditions into the User Guide.

MISC

* Moved info in MARSSsettings.R to .onLoad function.
* Added suppressWarnings() wrapper to KFAS call when R=0 in MARSSkfas since update to KFAS package produces warning messages when R=0.
* Typo in Eqn 124 of EMDerivation.pdf.  \beta should have been ^{-1}. Typo in Eqns 133 and 134.  vec parentheses should have been in front of R in second summation.  R in first line of eqn 133, was not referring to R (the var-cov matrix).  It should have had a new symbol.  Switched to T.  Eqn 134 was not R but this 'T'.
* Added safe to control list in man file MARSS.Rd   Left off accidentally.
* Small change to DLM chapter to clarify that rotation matrix only exists if Z has more than 2 columns.
* All subfunctions for a function moved into the main functions so they are hidden to the rest of the functions.
* The logLik function was using the logLik, samp.size and df attributes from the MLE object, but this is prone to creating errors.  The user may have changed the model structure or data in a MLE object and is trying to get the new logLik.  Changed to recompute the logLik.
* Removed use of stringr package; did not need
* Poor name choice. y.se was not the standardard error of ytT since sqrt(OtT) not sqrt(OtT-ytT^2) was being returned.  Changed name of y.se to ytT.se


MARSS 3.9 (CRAN 3-21-2014)
------------------------------------
none. resubmission due to missing file


MARSS 3.8 (CRAN 3-18-2014)
------------------------------------

ENHANCEMENTS

* Added check for fun.kf value in `checkMARSSinputs()`
* Added check to `print.marssMLE()` to make sure models are class marssMODEL.  
* Changed `summary.marssMODEL()` to return the list matrix instead of the marssMODEL passed in. Added tinitx to the returned (and printed) list.
* Removed `is.blockunconst()` and `is.blockequaltri()` functions.  Not really used or useful and were buggy.
* Much of this function code (assoc with identifying blocks) incorporated into a better is.validvarcov function to test for many more illegal constraints on a variance-covariance matrix.  This will catch most but not all illegal constraints on Q, R and V0.  It has a method argument, so method=BFGS can be passed in to check that all blocks are diagonal or unconstrained as needed by the chol transformation used in the `MARSSoption()` code to ensure varcov matrices stay postitive-definite.
* In `MARSS()`.  
   - Switched to use `MARSSkf()` to return kf (so use what user requested), but set Innov, Sigma, Kt etc with `MARSSkfss()`.
   - Added row names to states.se and y.se.
* In `MARSSkem()`. Removed adding of kf and Ey when trace>0.  This happens in MARSS().
* Changed `summary.marssMODEL()` to use marssMODEL attributes for par.names and model.dims, so it works on non-marss form marssMODEL objects.
* Added ability to handle time-varying var-cov matrices in     `MARSShessian()`
* Added check that Hessian CIs are only computed for models with diagonal var-cov matrices
* Added ability to deal with NAs in Hessian in `MARSShessian()`

DOCUMENTATION

* The OR and CA years for the harborSeal dataset were off.  Fixed and added references to the man file for harborSeal.  Removed the harborSealnomiss dataset as that is no longer used in the User Guide.
* Rewrote the Seal Population Structure chapter and MAR(p) chapter.
* Added info to `MARSSinfo()` to give the user some code to convert pre-3.5 marssMLE object to the 3.5+ form.

BUGS

* In `MARSSkfss()`. When Z was not square (num rows > num cols), OmgRVtt was not getting set.  OmgRVtt sets Vtt diagonals (and corresponding cols and row) to zero when R had 0s on the diagonal.
* In `MARSSkfas()`. Was returning $Innov and $Sigma using $v and $F, but as detailed in the KFS help page (KFAS package), the ones returned by KFS are not the same as the standard innovations and Sigma for multivariate data.  Now, `MARSSkfas()` returns a text message to use `MARSSkfss()` to get these.
* `residuals.marssMLE()` and `MARSSinnovationsboot()` were not running `MARSSkfss()` to get Innov, Kt, and Sigma when R was not diagonal.  Problem occurred after I changed `MARSSkfss()` to return text error instead of NULL for these.
* Bug introduced in 3.6 that printed no absol convergence when convergence=10.  Should have printed abstol convergence only.
* Bug in MARSSoptim (method=BFGS) that lead to only diagonal var-cov matrices when anything other than a diagonal var-cov matrix was selected.
* Same bug affected attempt to compute CIs for non-diagonal var-cov matrices with Hessian.
* Bug in MARSSoptim (method=BFGS) that allowed user to specify time-varying Q and R models, which code does not allow because cannot backsolve for par in that case.
* Bug in MARSSoptim (method=BFGS) that allowed Q, R, and V0 structures that can't be handled by the chol transformation in that code.  The transformation requires that Q, R, and V0 matrices be block unconstrained.  Blocks can be identical or unique or some identical and others unique but each must be unconstrained.  Note, in the context of a "block" matrix, a diagonal matrix is composed of n 1x1 blocks where n=nrows.  Thus by definition, a diagonal matrix (with shared or unshared elements on the diagonal) is always block unconstrained.  Dealt with with new is.validvarcov() function.
* Bug in convert.model.mat() when user used names like "2" or "1" and had fixed values of the same (e.g. 1,2).  This is because, inexplicably, R considers 1=="1" to be TRUE (and 2=="2", etc).  Replaced with sapply and identical() embedded within.
* There is a check in MARSSkfss() that any 0s on the diagonal of Vtt1 have a corresponding 0 on the diagonal of Q.  Was this line: Q0s=identical(which(diag.Q==0),which(diag.Vtt1==0)).  But that forced a more stringent requirement, that all 0s on diag of Q were identical to 0s on diag of vtt1 rather than that 0s on diag of Vtt1 had 0 on diag on Q, but not the converse.  Changed to Q0s=all(which(diag.Vtt1==0)%in%which(diag.Q==0)) so that requirement is one-way.
* X names were not getting applied to states in MARSS(); default X.names would be odd for non-design Z matrices. MARSS_marss() and MARSS_marxss().


MARSS 3.7 (CRAN 12-14-2013)
------------------------------------
Version 3.7 update was required due to new version of KFAS that changed its API.

ENHANCEMENTS

* Changed dependency to new version of KFAS.  Updated NAMESPACE to only import the 3 KFAS functions used.
* Created a versiontest.R file for comparing output from two different versions of MARSS.  It's in the doc directory.
* Exported the toLatex method for making a LaTex version (and pdf) of a marssMODEL object.

DOCUMENTATION

* Changes to the documentation. Mark made a few changes to the DFA chapter.  Got rid of Bluegreens in the example since that was mostly missing data.  Moved the info on the general MARSS equation from the Quick_Start guide to the chapter on algorithms.  Added a 'Tips and Tricks' section to the Quick_Start guide.
* Cryptomonas misspelled in data files
* The index file was out of date with old names for R scripts. Made some minor updates to the MARSS man file.

BUGS

* Fixed `allow.degen()` bug that would set elements to zero, leading to non positive definite matrices. Test if Q and R are diagonal.  If not, don't allow degens to be set since that is likely to lead to non-pos def matrices.  I could test if the row/col covariance are 0s but that would be costly.
* Fixed `loglog.conv.test()` bug that returned NAs when logLik > 720 due to exp(LL) call.  Changed to exp(LL-mean(LL))


MARSS 3.6  (CRAN 11-26-2013)
------------------------------------
Version 3.6 update is mainly concerned with speeding up `MARSS()` for problems with large number of time series (n > 100) and where many R elements are being estimated (e.g. R="diagonal and unequal").  This comes up in dynamic factor analyses often.  The changes also improve speed for small R problems by about 25%, but speed increase is 10 fold for problems with R matrices that are 100x100 with 100 estimated R elements.

ENHANCEMENTS

* Changed `MARSS.marxss()` to speed up conversion of "unconstrained" shortcut to a matrix.  Only matters if m or n is big.
* Sped up `convert.model.matrix()`.  Old version was always using slow code to deal with * and + in character matrix.  This made formation of the free and fixed matrices very, very slow when the matrix got big (100x100, say).
* Added silent==2 which gives verbose output of progress
* Changed `is.design()` to not use near equality for test of element==0.  This may break `MARSS()` since R sometimes doesn't maintain "zeroness".
* Removed many inefficiencies in `MARSSkem()` code for working with large matrices.  Replaced all crossproducts with `crossprod()` and `tcrossprod()` which are significantly faster for large matrices.  This increases speed 2-10 fold when working with larger matrices.  Largest speed increases are when R is not diagonal and equal.
* Hard coded a fast diagonal test into `MARSSkfss()` instead of using the very slow `is.diagonal()` function (which is really meant for list matrices)
* Added set.degen to `degen.test()` function so that it sets a flag to TRUE if any var-cov diagonals set to 0.  If so, do the updates otherwise skip.
* Improved speed of `parmat()` by testing if d and f matrices are not time-varying.  In which case, don't subset the array, but rather rest the "dim" attribute.  Much, much faster for big d and f matrices.
* Improved `sub3D()` to make it a bit faster by using x[,,t] when both nrow and ncol are >1
* Improved `vec()` to make it 3x faster by setting dim attr instead of using matrix() when matrix is 2D

DOCUMENTATION

* `lakeWAplankton` datasets were saved as data.frame.  Changed to matrix.
* Created R files for all the 'application' chapters in the user guide.

BUGS

* Fixed bug in building A matrix for A="scaling" which would throw warning if zero columns were in Z.  No error but just unnecessary warnings.
* `MARSS()` didn't print out marssMLE object when convergence=12 (maxit set below min for conv test).
* If there were duplicated rownames in the data, R and U would use those to set shared values.  This was a bug.  Added a test for duplicated rownames, and add "-1", "-2" etc to duplicated name to distinguish them.

MARSS 3.5
------------------------------------
Version 3.5 is mainly concerned with formalizing the internal structure of model objects.  marssMODEL objects have been formalized with attributes.  A form definition along with associated form functions have been defined.  This won't be noticeable to users but makes writing functions that use marssMODEL objects easier and more versatile.

ENHANCEMENTS

* changed `MARSS.dfa()` to allow B and Q setting to "diagonal and equal" or "diagonal and unequal"
* Fixed the printing of model structure so it shows the form that the user called rather than base form
* Added basic predict function (note, not exported to users in 3.5.  accessible to users via `MARSS:::predict.marssMLE()`. Further development will be done before exporting to users.)
* Removed `MARSSvectorizeparam()` and `MARSSapplynames()` from the exported list.  The former has been replaced by coef(marssMLEObj, type="vector").  The latter is an internal utility function.
* Changed `MARSSkfas()` to return Innov and Sigma when R is diagonal.  When R is not diagonal, the user is directed to use `MARSSkfss()` since `MARSSkfas()` and `MARSSkfss()` do not agree when R is not diagonal (and I think the error is in KFAS as the Sigma looks off when R not diagonal).
* Changed `MARSShessian()` to use a Cholesky transformation on any variances so that the variance covariance matrices stay positive definite
* Change above required update to `MARSSparamCIs()`
* miss.value is now deprecated.  The user is instructed to replace missing values with NAs before passing data to `MARSS()`.
* Created an global environment (pkg_globals) specific to the package environment, so that all functions have access to these package-specific globals.  This is assigned in a new onLoad() function.
* Added check in `MARSSkfas()` for version of KFAS.  API changes in KFAS 1.0.0, and a line of code was added to use the correct API if KFAS 0.9.11 versus 1.0.0 is installed.  MARSS will work with both versions of KFAS.
* Condensed the errors print-out to 10 errors; added more error info to `MARSSinfo()`
* Restructured the NAMESPACE and DESCRIPTION files to better control imports and dependencies so that user cannot break the package by detaching needed libraries or redefining needed base and stats functions.

DOCUMENTATION

* Added better message reports when model list elements are not allowed
* Updated the help files
* Added more items to `MARSSinfo()`
* The original Lake Washington dataset was added as `lakeWAplanktonRaw`.  Month^2 dropped and month not z-scored.  Original raw data for the counts.
* Added dynamic linear model case study
* Revamped and extended the covariates case study

BUGS

* `MARSSboot()` was out of date with newest version of `MARSShessian()`'s returned arguments.
* `is.blockunconst()` bug made it break on certain diagonal list matrices

MARSS 3.3 and 3.4  (CRAN 1-16-2013)
------------------------------------
This version update is mainly concerned with adding generic functions (coef, residuals, predict), hooking back up KFAS package filters into MARSS functions, and customizing print functions for different model forms.

ENHANCEMENTS

* Linked the KFAS package to MARSS
* `MARSSkfas()` was changed to work with the new KFAS version released July 2012.  This led to a 10-20 fold decrease in computation time for method="BFGS" and 2 fold for method="kem".
* `MARSSkf()` changed to `MARSSkfss()`; `MARSSkf()` is now a utility function that picks `MARSSkfas()` or `MARSSkfss()` based on MLEobj$fun.kf
* Added a lag-one covariance smoother to the output of `MARSSkfas()` using the algorithm given on page 321 in Shumway and Stoffer (2000), Time Series Analysis and Its Applications (note the 2000 edition not 2006).  The algorithm is given in the User Guide also.  The EM Algorithm requires the lag-one covariance smoother but this is not one of the outputs of the `KFS()` function in the KFAS package.
* Changed the Kalman filter output in `MARSSkfas()` to be strictly 0 when it is supposed to be (when R has 0 on the diagonals); `KFS()` output has ~0 not actually 0.
* Changed the print function for marssMLE objects so that printing can be customized for model form.
* Added `coef()` method for marssMLE objects.  Added $coef to marssMLE object.
* Changed `parmat()` to be hidden (not exported).  Instead its functionality is through the standard R function for this purpose, `coef()`.
* Added `residuals()` method for marssMLE objects by changing `MARSSresids()` to `residuals.marssMLE()`. 
* Added `predict()` method for marssMLE objects.
* Added $call, $alt.forms to marssMLE objects so that printing (and coef and other functions for model objects) can be customized to form.
* Added the standard error of the missing y values.  MLEobj$y.se.  Edited .Rd file to reflect changes.
* Added utility function all.equal.vector to test for equality in vectors and matrices
* Changed fixed.free.to.formula to allow 3D matrices; returns array if fixed/free indicated time-varying parameter
* Added `toLatex.marssMODEL()` function to create latex and pdf output of models
* Various changes to NAMESPACE in conjunction with the above changes.

BUGS

* `MARSS.dfa()` did not allow user to pass in Z as matrix.
* `parmat()`. When t was a vector, `parmat()` only returned the value at max(t).
* `MARSSkemcheck()` crashed when the test "when u^{0} or xi^{0} are estimated, B adjacency matrix must be time invariant" was started.
* `MARSS_marxss()` threw an error when Z was passed in as a matrix and A="scaling"
* `describe.marss()` bug caused diagonal matrices with 1 estimated value and fixed values to not be identified as diagonal

MARSS 3.2  (CRAN 08-28-2012)
------------------------------------
Version 3.2 is a minor update to the documentation

DOCUMENTATION

* Some edits to the case studies and the User Guide to fix typos and stuff noted in the August workshop
* Added data to the Isle Royale dataset including covariates (temperature and precipitation)
* Added isleRoyal.Rd man file for Isle Royale data and covariates.
* Fixed misspelling in DESCRIPTION

BUGS

* Fixed bug that prevented `MCInit()` from working.

OTHER

* Moved .Rinstignore to the top-level


MARSS 3.0  (CRAN 07-10-2012)
------------------------------------
Version 3.0 is a major update and clean-up. Besides the clean-up, the changes were to allow time-varying parameters and a way for the user to specify linear constaints using an eqn like `a+2*b` in the parameter matrix.

The changes are extensive but are internal and should be largely invisible to users of MARSS 2.X.  The `MARSS()` 3.0 call is backwards compatible to 2.9 except that kf.x0 changed to tinitx and moved from control list to model list.   Use of KFAS remains disabled until I can update to the new version of KFAS.  This slows down method="BFGS", but does not affect method="kem".

INTERNAL CHANGES

* Meaning of fixed and free changed. fixed is the f matrix in the EM derivation and free is the D matrix.  Originally, I used a fixed/free pair with NAs.  The new form closely follows the derivation and leads to more unified code.  This required changes in most files to deal with new meaning of fixed and free.
* Added element X.names to marssm model object
* Added model list to model object so that what the user passed in with the `MARSS()` call is retained (for reference).
* Allow user to spec linear constraints. The main point of 3.0 is so that users can spec intercept + beta_1 p_1 + beta_2 p_2 ... constraints.  Users do this with the usual list matrices using something like "theta+phi" or "2+2\*theta+phi".  I added a function to interpret basic math like that with + and \*.
* Changed args for `MARSSkf()` and `MARSShatyt()` functions. Takes MLEobj now.
* Added `parmat()` function. This returns a parameter matrix when given a MLEobj.
* Changed par element of MLEobj. It is a vector of only the estimated elements now.  This required changes to user manual where I show specific parameters. If par is fixed, par element is matrix(0,0,1).
* Changed kf.x0 and t.x0 to tinitx. This standardized the naming a bit and the name now stands for "t at initial x" which will hopefully be easier to remember.  Removed "x00" and "x10" at least where user will see it.  They are still in the internal code. tinitx is passed in in the model list.
* Changed the way the marssm objects are created. Before the code locked one into a MAR-1 state-space form.  However many different types of time series models can be rewritten in MAR1SS form.  This rewriting is onerous for users and I don't want them to have to do that.  Also I wanted to make it easier to write functions to write different time series models in MAR-1 SS form.  Now `MARSS()` looks for a function called MARSS."form", where form is something like "mar1ss" or "dlm".  This function takes the MARSS inputs (all of them) and transforms the input into a marssm object ready for the fitting functions as the function writer wishes.  All that the function has to do is to return a valid marssm object from the model element of the `MARSS()` call.  This allows me (or anyone else) to use whatever parameter names they want in the model element.  This way the user can use familiar names for parameters can set some parameters to specific values (like 0).  Or the user could do something totally different with the model element and just have it be a text string like model="my.ts.model.1" or model="my.ts.model.2".   The only constraint is that the function output a proper marssm object and that the control, inits, and MCbounds arguments for MARSS are properly specified.
* Removed popWrap.r, checkPopWrap.r, MARSSoptions.r.  Became obsolete with above changes
* Added `checkMARSSInputs()` and `checkModelList()`. These replaced the functionality of popWrap.r and checkPopWrap.r
* Added `MARSS.marxss()`.  This is the first `MARSS.form()` function.  This is a standardized format for so that I can add other forms easily.
* Changed `MARSSkf()` so that K (Kalman gain) is 0 when tinitx=1 and V0=0. Changed `MARSSkf()` to allow some of diagonals of V0 to be 0 and others non zero.  Got rid of many of the OMGs.  Added `pcholinv()` function to diaghelpers.r which deals with matrices with 0s on diagonals.  This streamlined the filter code.
* Rewrote many sections of `MARSSkem()` to allow time-varying parameters. Made changes to `MARSSkf()`, `MARSSkfas()` and `MARSSsimulate()` to allow time-varying parameters.  See EMDerivations.pdf
* Added fun argument to `MARSShessian()` and `MARSSparamCIs()` to allow one to specify the function used to compute the log-likelihood.
* Added row and col names to Hessian in `MARSShessian()`
* Moved diffuse from control element to model element of MLEobj since it is part of the model specification.  Required changes to `MARSSsettings()`, `MARSS.marxss()`, `is.marssm()`, `is.marssMLE()`.
* Changed `MARSSkem()` and `MARSShatyt()` to allow some diag.V0=0 and others not 0, so user can mix stochastic and fixed initial x states.
* Rewrote x0 and U update sections in `MARSSkem()`.  Removed the OMGs from `MARSSkem()` since no longer needed given the new `pcholinv()` function.

DOCUMENTATION

* Totally revamped the EMDerivations.pdf to allow time-varying models.
* Rewrote (again) the section on degenerate models in EMDerivation.pdf to allow B structures that imply both total deterministic X and indirectly stochastic x.  The latter is required to allow one to rewrite a MAR-p model as a MAR-1 model.  Time-varying params meant that the matrix geometric function no longer could be used, but I found a simplier recursion.  Improved the presentation so only 1 x0 and U update equation is given rather than 5 special cases.  

OTHER CHANGES

* describe_marssm, rewrote
* MARSSmcinit, changed how did draws
* MARSSparamCIs, rewrote, changed how I store se, upCIs, etc.  now as vector like paramvec
* MARSSvectorizeparam, rewrote
print and summary functions, updated
* MARSSinits, rewrote, returns new form of parlist
* MARSSkem, changes to R, Q, x0 & U update per new degenerate model update eqns, added p, removed fixed and replaced with f, removed free and replaced with d
* is_marssm, added X.names to model
* as_marssm, removed and replaced with MARSS.marxss
* Removed fixed and free from allowable MARSS() input (affected MARSSsettings, PopWrap, popwrapcheck)
* MARSSLLprofile, removed for now, not sure it works
* MARSSoptions, removed, obsolete
* MARSScheckdims, removed, not used
* MARSScheckpar, removed, not used
* popWrap and checkPopWrap, removed, functionality replaced with checkMARSSInputs and CheckModelList
* diaghelpers.r, added parmat, pcholinv, pinv, few other functions

BUGS

* Bug in `MARSSkem()` that meant that the maxit-1 kf and logLik were returned when algorithm stopped due to hitting maxit.  Par was correct.


MARSS 2.9 (2012-03-20)
------------------------------------
Version 2.9 was a temporary update to deal with a major change to the API of the KFAS package. Needed to disable use of `MARSSkfas()` until that function was rewritten.

ENHANCEMENTS

* Fixed `MARSSboot()` so that MLE objects with method=BFGS can be used; changed the param.gen argument to take "MLE" and "hessian" instead of "KalmanEM" and "hessian".  Updated MARSSboot.R and MARSSboot.Rd.
* Temporarily disabled calls to `MARSSkfas()` until MARSS can be made compatible with new version of KFAS package.  Removed importFrom(KFAS, kf) and importFrom(KFAS, ks) from the NAMESPACE.  Removed MARSSkfas from the export list in NAMESPACE.  Removed KFAS in the depends line of DESCRIPTION.

DOCUMENTATION

* Updated the DFA example in the manual.
* Changed the column headings in the L WA plankton dataset slightly to have uniform capitalization.

BUGS

* Fixed `MARSSaic()` and `MARSSparamCIs()` so that `MARSSboot()` call uses param.gen="MLE".  This fixes the bug that stopped MLE objects from BFGS calls to fail.


MARSS 2.8 (2012-01-23)
------------------------------------
Version 2.8 improved default initial conditions functions and fixed bugs in the Shumway and Stoffer Kalman filter/smoother function.

* Added NEWS file, .Rinstignore in inst\doc
* Added example of lag-p model to the manual.
* Fixed bug in `MARSSkf()` when R=0, kf.x0=x10, and V0=0. The algorithm was not setting x(1) via y(1) in this special case.
* In `MARSSinits()`, got rid of the linear regression to get inits for x0; using instead solution of pi from y(1)=Z\*(D\*pi+f)+A; This stops MARSS from complaining about no inits when Z is not a design matrix.  NOTE NB: This means the default initial x0 are different for 2.7 and 2.8, which leads to slightly different answers for `MARSS(dat)` in 2.7 and 2.8. The answers are not really different, just they started with slightly different initial values so have slightly different values when the algorithm reaches its convergence limit.
* Removed dependency on time package. The progressBar function was moved into MARSS since the time package is no longer maintained.
* Changed `MARSSkemcheck()` to allow lag-p models. I worked on the derivation of the degenerate models (with 0 on diag of Q) to better define the needed constraints on B.0 and B.plus sub matrices.  This led to changes in MARSSkemcheck.r so that lag-p models written as MARSS model are now allowed.  There are still problems though in x0 estimation in the EM algorithm when there are zeros on R and B diagonals, so best to method=``BFGS'' until I redo the degenerate EM algorithm.
* Added option to force use of `MARSSkf()` function instead of MARSSkfas. If kf.x0="x10", default was to use `MARSSkfas()` function which is much faster, but it doesn't like 0s on B diagonal if V0 is 0.  So I added the option to force use of slower `MARSSkf()` function using method="BFGSkf". Reguired adding stuff to MARSSsettings.r and MARSSoptim.r.  This is mainly for debugging since `MARSSoptim()` will now check if optim failed and try using `MARSSkf()` if `MARSSkfas()` was used.  Added line to output that says which function used for likelihood calculation; again for debugging.
 * Edited `MARSSmcinit()` to improve random B generation. There is nothing to guarantee that random Bs in mcinit routine will be within the unit circle, however it is probably a good idea if they are.   Default bounds for B changed to -1,1 and random B matrix rescaled by dividing by max(abs(Re(eigen(B)))/runif(1) to get the max abs eigenvalue between 0 and 1.  This works unless the user has fixed some B values to non-zero values.  This required change to is\_marssMLE.r also to remove constraint that B bounds be greater than 0.
* Edited `MARSSmcinit()` to allow fixed and shared values in random Qs and Rs. The random Wishart draw is rescaled based on the fixed and shared structure in R or Q.  As part of this, I cleaned up how fixed and shared values are specified in the random draws of other parameters.  This change doesn't change the end effect, but the code is cleaner.

BUGS

* `MARSSoptim()` did not allow unconstrained Q or R. The problem had to do with temporarily resetting the upper triangle of tmp fixed matrices to 0 when using tmp.par as chol matrix.
* Error in `MARSSkf()` when there were 0s on diagonal of Q. The algorithm only worked if B was diagonal.  Fix required changes to Kalman smoother bit of `MARSSkf()`. I rewrote the pertinent section in EMDerivation.pdf.

DOCUMENTATION

* Cleaned up degenerate derivation in EMDerivation.pdf
* Added warning in the covariate section. The error-free covariate section in the manual did not clarify that the log-likelihood of the covariates with the dummy state model would be included in the MARSS output.  MARSS version 2.9 will allow error-free covariates in a more standard manner.


MARSS 2.7 and 2.6 (2011-10-21)
------------------------------------
Versions 2.7 and 2.6 focused on misc. bugs.

* Added sections on covariates and lag-p model to the user guide.
* MCInit was not working for non-diagonal R and Q. I replaced the function for randomly drawing matrices with a random draw from a Wishart distribution.
* m not getting assigned in `MARSSPopWrap()`.  Some of the allowable cases for Z and m were missing.
* Added more info re R or Q not positive-definite in error messages. If the user specifies an illegal variance covariance structure from a general estimation perspective (nothing to do with MARSS), they can get the "not positive-definite" error.  Added some text in Troubleshooting section to help if they get this error.
* Fixed `MARSSsimulate()` bug. `MARSSsimulate()` was broken for multivariate simulation since I forgot that rmvnorm returns a 1 x p matrix even if the mean is p x 1.  Wrapped the rmvnorm call in a array() to fix the dim setting.
* Error in x0 update when R=0 and x0 fixed. If x_1 has fixed elements, estimates should not used for those elements.  Code was missing some d$x0 bits.  This means that the user can fix x_1 when R=0 to a value not equal to the corresponding y_1 value.  This would mean an illogical model so a check was added to stop and give warning if that happens.


MARSS 2.5
------------------------------------
Version 2.5 focused on switching model specification to use list matrices.

* Factor option for all but Z removed. Same functionality is now provided via list matrices
* Removed fixed/free args from `MARSS()`. Same functionality is provided via list matrices
* Constraint arg changed to model in `MARSS()`. Just the name of the argument was changed to be more intuitive
* Rewrote user guide to reflect above changes
* Added case studies to user guide on dynamic factor analysis and species interactions with covariates


MARSS 2.2
------------------------------------
Version 2.2 focused on incorporating the KFAS Kalman filter/smoother functions which are faster and more stable.

* Added diffuse priors for method="BFGS" and kf.x0="x10"
* Incorporated KFAS package. Their Kalman filter is faster but only for x10.  Added `MARSSkfas()` function.
* Changed Q/R estimation in optim to allow off-diagonal terms.
* Added V0 estimation option.  This works like other parameters now
* LL calc when R=0 fixed. LL calc in `MARSSkfas()` to deal with 0s on diag of Ft[,,1] so can do R=0
* Replaced `show.doc()` with `RShowDoc()` (base)
* Default miss.value changed NA where NA is as.numeric(NA) rather than logical.


MARSS 2.0
------------------------------------
Version 2.0 implements changes to allow B and Z estimation and more element sharing in Q and R matrices.

ENHANCEMENTS

* `MARSSkem()` algorithm changed to allow B and Z estimation.
* `MARSSkem()` algorithm changed to allow constrained B and Z estimation.  This was the second main objective of MARSS 2.0.  This allows you to have fixed values or shared values in your B or Z matrices.
* Allow more types of element sharing in the Q and R estimation. In MARSS 1.1, you were limited to diagonal, equal var-cov, and unconstrained.  Now various types of block-diagonal matrices are allowed.
* Allow some Q or R variances to be set to 0. This allows partially deterministic systems (Q=0) and systems with no observation error (R=0)
* Fixed the V0=0 case. I was using a work-around to do the fixed x at t=0 case (V0=0).  I derived the solution and added this to `MARSSkem()`.  There is no iter.V0 control element anymore.
* Changed logLik conv test. I was doing the log-log test against logLik instead of log(logLik).  I think the test works better using the log of the log-likelihood.
* Detect degeneracy and set Q or R element to zero. Now instead of the variance walking to log(negative infinity) in an infinite number of iterations, the algorithm detects that a variance is going to zero and tries setting it to zero.
* `MARSSkem()` changed to a more general way to deal with missing values. This is described in the EMDerivation.pdf.  It doesn't affect the user, but allows the code to be expanded to more types of models much more easily.
* Changed to using list matrices to describe models. Now you can essentially write the way your model looks on paper (in matrix form) as a list matrix in R and it will run.  No more fixed and free matrices---at least from the user's perspective.
* Added some code optimization. I cleaned up some of the things that really slowed down 1.1.  2.0 is now about as fast as 1.0 was.

DOCUMENTATION

* Big revamp of EMDerivations.pdf. I cleaned up my derivation a lot.  I'm especially happy with the sections on dealing missing values part of the derivation.  It's much more elegant and logical now.  The sections on degenerate matrices are cluttered and the notation is painful, but I will leave them be for awhile.

BUGS

* Bug in miss.value=NA. When miss.value=NA, class for NA was logical.  Needed to be numeric.


MARSS 1.1
------------------------------------

* Fixed formatting issues with error messages. 
* Allow NA and NaN to be used for miss.value
* Fixed bug in `MARSSmcinit()`. MCMC init function would crash for anything except the default model.
* Fixed ungraceful exiting when minit > maxit
* Fixed ungraceful exiting when method=BFGS threw error
* Added more info to ?MARSS and help(``MARSS-package''). Changed MARSS.Rd and MARSS-package to have reference to user guide, index, and MARSS-package help page.
* Changed convergence test. In the convergence diagnostics test, we check that the slope of logLik vs (log iteration num) is close to zero.  This is a standard convergence test.  But Shumway and Stoffers code uses a delta logLik test which checks that the logLik.new-logLik.old is less than some absolute (user specified) tolerance.  This turns out to be a bad convergence test because the log-log plot (described above) can still have a fairly clear slope.  I switched over to using the log-log test as the default test, but I allow the user to specify a abstol (delta logLik) if they want that instead.  This change slows down model fitting considerably but model fits that are actually converged.\
* Fixed `is.design()` function. A design matrix must have more or equal rows than columns.
* R was changing dims on some matrices in `MARSSkf()`. R has a flaw in terms of how it behaves when you subscript a matrix and the new matrix has a dimension length of 1 for one (or more dimensions).  For example, if a=array(0,dim=c(1,2,4)), then a[,,1] is no longer a matrix but instead is a vector and dim(a[,,1]) is NULL.  This can cause all sorts of mysterious bugs.  Sometimes adding drop=FALSE will prevent this unpleasant behavior.  If b=matrix(0,2,2), dim(b[,1,drop=FALSE]) is c(2,1) while dim(b[,1]) is NULL.  drop=FALSE works great with 2-dimensional matrices, but with 3-dimensional matrices it doesn't work.  If a=array(0,dim=c(1,2,4)), dim(a[,,1,drop=FALSE]) is c(1,2,1) instead of c(1,2) which is what you want if a[,,1] is what is going to appear in some matrix operation. This problem came up in the Kt[,,t] %\*% innov[,t] line in MARSSkf.  Normally Kt[,,t] is square and a square matrix or a scalar is returned, but if Kt[,,t] happened to be something like dim=c(1,3,20) then Kt[,,t] returned a VECTOR of length 3.  In this case, Kt[, , t] %\*% innov[, t] crashed the code.  I had to use a kluge to force R to keep the dimensions after subscripting. This bug only occurred in models where Z is not a design matrix.
* Fixed formatting issues in summary(marssm object). The naming of elements in the model matrices did not match `summary(marssMLE object)`.
* Added function `MARSSoptions()`. This allows you to change the defaults for the `MARSS()` function.  See `?MARSSoptions`.
* Added function `MARSSLLprofile()`. This allows you to plot some basic log-likelihood profiles.  See ?MARSSLLprofile.