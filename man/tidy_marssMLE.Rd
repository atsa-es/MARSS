\name{tidy.marssMLE}
\alias{tidy.marssMLE}
\title{Return estimated parameters, expected value of X(t) and Y(t) with summary information}
\description{

\code{tidy.marssMLE} returns summary information about the model parameters and estimated state and observation processes. In all cases, all the data are used, thus conditioning is on the data from 1 to T.

For parameters, \code{tidy.marssMLE} returns their estimates and their confidence intervals. For states (X) and observations (Y), it returns the expected values (mean value) and intervals. If you want to analyze your model residuals or are doing a crossvalidation (leave-one-out or k-fold), you need the residuals intervals. These are the intervals for model residuals (data - model fitted value). Use \code{\link{augment.marssMLE}} for the model residuals and their intervals.

The tidy function is designed to work with the \code{broom} package and you will need to load that package if you want to call \code{tidy(fit)} instead of \code{tidy.marssMLE(fit)}. 
}
\usage{
  tidy.marssMLE(x, type = c("parameters", "xtT", "fitted.ytT", "ytT"),
                conf.int = TRUE,
                conf.level = 0.95,
                form=attr(x[["model"]], "form")[1], ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type}{What you want estimates and intervals for. Parameters, smoothed states (xtT), the fitted y (Z xtT + A), or observations conditioned on all the data (ytT). If you want intervals for new data sets, use fitted.ytT. If you are getting estimates for missing data, use ytT. See details.}
  
  \item{conf.int}{ Whether to compute confidence and prediction intervals on the estimates. }
  
  \item{conf.level}{ Confidence level. alpha=1-conf.level}

  \item{form}{If you want the tidy function to use a different form than that specified in attr(x$model, "form").  Useful if you have a DFA model that you manually set up, which does not have the form attribute set. Normally just ignore and let the function use the "form" set in the attributes.}
  
  \item{...}{Optional arguments.  If \code{conf.int=TRUE}, then arguments to specify how CIs are computed can be passed in. See details and \code{\link{MARSSparamCIs}}. If form="dfa", \code{rotate=TRUE} can be passed in to rotate the trends (only trends not Z matrix).}

}

\value{
A data.frame with estimates, sample standard errors, and confidence (or prediction) intervals.
}

\details{

Below, X and Y refers to the random variable and x and y refer to a specific realization from this random variable.

\strong{type="parameters"}  

If \code{type="parameters"}, this returns a data.frame with the estimated parameters  of a MARSS model with, optionally, standard errors and confidence intervals. This assembles information available via the \code{\link{print.marssMLE}} and \code{\link{coef.marssMLE}} functions into a data.frame that summarizes the estimates.  If conf.int=TRUE, \code{\link{MARSSparamCIs}} will be run to add confidence intervals to the model object if these are not already added.  The default CIs are calculated using a analytically computed Hessian matrix.  This can be changed by passing in optional arguments for \code{\link{MARSSparamCIs}}. 
  
\strong{type="xtT"}

\code{tidy.marssMLE} returns the confidence and prediction intervals of the state at time t conditioned on all the data and using the estimated model parameters as true values. The prediction intervals (and .sd.x) are the standard intervals that are shown for the estimated states in state-space models. For example see, Shumway and Stoffer (2000), edition 4, Figure 6.4. As such, this is probably what you are looking for if you want to put intervals on the estimated states (the x). However, these intervals do not include parameter uncertainty. If you want state residiuals (for residuals analysis), use \code{\link{residuals.marssMLE}} or \code{\link{augment.marssMLE}}.

\strong{Quantiles} The state \eqn{\mathbf{X}_t} in a MARSS model has a conditional multivariate normal distribution, that can be computed from the model parameters and data. In Holmes (2012, Eqn. 11) notation, its expected value conditioned on all the observed data (1:T) and the model parameters \eqn{\Theta} is \eqn{\tilde{\mathbf{x}}_t}{\tilde{x}_t}. In \code{\link{MARSSkf}}, this is \code{xtT[,t]}. The variance of \eqn{\mathbf{X}_t}{X_t} conditioned on the observed data and \eqn{\Theta} is \eqn{\tilde{\mathbf{V}}_t}{\tilde{V}_t} (\code{VtT[,,t]}). Note that VtT[,,t] != B VtT[,,t-1] t(B) + Q, which you might think by looking at the MARSS equation for x. That is because the variance of W(t) conditioned on the data (past, current and FUTURE) is != Q (Q is the unconditional variance).

\eqn{\tilde{\mathbf{x}}_t}{\tilde{x}_t} (xtT) is an estimate of \eqn{\mathbf{x}_t}{x_t} (the true value), and the standard error of that estimate is given by \eqn{\tilde{\mathbf{V}}_t}{\tilde{V}_t} (\code{VtT[,,t]}). Let \code{se.xt} denote the sqrt of the diagonal of \code{VtT}. The equation for the \eqn{\alpha/2} confidence interval is (\code{qnorm(alpha/2)*se.xt + xtT}). \eqn{\mathbf{x}_t}{x_t} is multivariate and this interval is for one of the \eqn{x}'s in isolation. You could compute the m-dimensional confidence region for the multivariate \eqn{\mathbf{x}_t}{x_t}, also, but \code{tidy.marssMLE} returns the univariate confidence intervals. 

The variance \code{VtT} gives information on the uncertainty of the true location of \eqn{\mathbf{x}_t} conditioned on the observed data. As more data are collected (or added to the analysis), this variance will shrink since the data, especially data at time t, increases the information about the locations of \eqn{\mathbf{x}_t}. This does not affect the estimation of the model parameters, those are fixed, but rather our information about the states at time t.

If you have a DFA model (form='dfa'), you can pass in \code{rotate=TRUE} to return the rotated trends.  If you want the rotated loadings, you will need to compute those yourself:
\preformatted{
dfa <- MARSS(t(harborSealWA[,-1]), model=list(m=2), form="dfa")
Z.est <- coef(dfa, type="matrix")$Z
H.inv <- varimax(coef(dfa, type="matrix")$Z)$rotmat
Z.rot <- Z.est \%*\% H.inv
}

\strong{Intervals for the observation process}

For observation process, the expected values and intervals are shown for either new data (\code{type="fitted.ytT"}) or the observed data set  (\code{type="ytT"}). Details on these are below after this discussion of intervals for the observation process

The types of intervals you want for data (Y part of the MARSS equation) depends on what you are trying to do.
\itemize{
\item{Get the model predictions of the expected value of new Y or some underlying mean Y}{Use \code{type="fitted.ytT"}. This returns the fitted values (model predictions = Z x(t)+A) for \eqn{Y_t} conditioned on all the data. Confidence intervals and prediction intervals are returned. The former is the interval for the mean of new data and the latter is the interval for new data (not the mean but data themselves).}
\item{Get the distribution of new data at time t that would be generated by the process}{Same as above.}
\item{Compare your data to model predictions}{In this case, you want the distribution of the model residuals for the data. Use \code{\link{augment.marssMLE}} with \code{type="observations"}. You want the standard errors for the observed data minus the fitted values which is what  \code{\link{augment.marssMLE}} gives.  }
\item{Get estimates and variance of missing data in your data set}{Use \code{type="ytT"}. The observed data will have an expected value equal to the observed data and variance of 0, while the missing data will have an expected value and variance conditioned on all the observed data. Note, if R is diagonal then the missing data values (and intervals) will be the same as for \code{type="fitted.ytT"} but if R is non-diagonal and some y at time t are missing and some are not, then the expected values will be very different.}
\item{Do a leave-one-out cross-validation}{In this case, you want the distribution of the model residuals for those left-out values. Use \code{augment.marssMLE} with \code{type="observations"}. You want the standard errors for the left-out data minus the fitted values which is what  \code{\link{augment.marssMLE}} gives.  }
\item{One-step-ahead predictions}{Use \code{fitted.marssMLE} with \code{conditioning="t-1"}. Confidence (mean prediction) and prediction intervals (new data) are returned.  }
\item{Y prediction conditioned on data up to t-1}{Same as one-step-ahead.}
}

\strong{type="fitted.ytT"}

For \code{type="fitted.ytT"}, \code{tidy.marssMLE} returns the analogous information for the Y part of the MARSS equation for an I.I.D. NEW DATA SET  \eqn{y'}. The expected value and variance of \eqn{y'} is conditioned on the data you did observe \eqn{y}. It is important to note that \eqn{y'} is independent and identical (meaning i.i.d. in a statistical sense) to \eqn{y} except it has no missing values. Do not plot your observed data on these intervals. You need residuals intervals in that case. See \code{\link{augment.marssMLE}} for those. 

The expected value of a new data set \eqn{\mathbf{Y'}_t}{Y'_t} conditioned on the observed data \eqn{\mathbf{Y}=\mathbf{y(1:T)}}{Y=y(1:T)}, is \eqn{Enewy = Z_t \hat{x}_t + D_t d_t + a_t}{Enewy = Z_t hatx_t + D_t d_t + a_t}, where \eqn{\hat{x}_t}{hatx_t} is the expected value of \eqn{\mathbf{X}_t}{X_t} conditioned on the data up to T. The variance of a new data set \eqn{\mathbf{Y'}_t}{Y'_t} conditioned on the observed data 1:T, is \eqn{var.newy = R_t + Z_t \hat{V}_t Z_t^\top}{var.newy = R_t + Z_t hatV_t t(Z_t)}, where \eqn{\hat{V}_t}{hatV_t} is the variance of \eqn{\mathbf{X}_t}{X_t} conditioned on the data up to T. The variance of the expected value of the new data set \eqn{\mathbf{Y'}_t}{Y'_t} is \eqn{var.Enewy = Z_t \hat{V}_t Z_t^\top}{var.Enewy = Z_t hatV_t t(Z_t)}.

We compute the prediction interval for y', an interval that will cover the new data for alpha/2 percent of new data sets.  The equation for the \eqn{\alpha/2}{alpha/2} confidence interval is (\code{qnorm(alpha/2)*sd.newy + Enewy}), where \code{sd.newy} is the square root of the diagonal of \eqn{var.newy}. The confidence interval for the expected value of y' is \code{qnorm(alpha/2)*se.Enewy + Enewy}, where \code{se.Enewy} is the square root of the diagonal of \eqn{var.Enewy}.

\strong{type="ytT" for missing data estimation}

This returns the expected value and variance of \eqn{\mathbf{Y}_t}{Y_t} (the data set you DID observe) conditioned on \eqn{\mathbf{Y}_t=y_t}{Y_t=y_t}. If you have no missing data, this just returns your data set. But you have missing data, this what you want in order to estimate the values of missing data in your data set. The expected value of \eqn{\mathbf{Y}_t|\mathbf{Y}=\mathbf{y}(1:T)}{Y|Y=y(1:T)} is in \code{ytT} in \code{\link{MARSShatyt}} output and the variance is \code{OtT-tcrossprod(ytT)} from the \code{\link{MARSShatyt}} output. 

The intervals reported by \code{tidy.marssMLE} for the missing values takes into account all the information in the data, specifically the correlation with other data at time t if R is not diagonal. Do not use \code{type="fitted.ytT"} for interpolating missing data as those are for entirely new data sets and thus will ignore relevant information if \eqn{\mathbf{y}_t}{y_t} is multivariate, not all \eqn{\mathbf{y}_t}{y_t} are missing, and the R matrix is not diagonal.

The standard error and confidence interval for the expected value of the missing data along with the standard deviation and prediction interval for the missing data are reported. The former uses the variance of \eqn{E[Y(t)]} conditioned on the data while the latter uses variance of \eqn{Y(t)} conditioned on the data. \code{\link{MARSShatyt}} returns these variances and expected values. See Holmes (2012) for a discussion of the derivation of expectation and variance of Y(t) conditioned on the observed data (in the section 'Computing the expectations in the update equations').

\strong{Parameter uncertainty} Currently the intervals calculations for the states and observations use the point estimates of the model parameters and thus solve the intervals for the 'known' parameters case.

}
\references{ 
R. H. Shumway and D. S. Stoffer (2000).  Time series analysis and its applications. Edition 4. Springer-Verlag, New York.

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]
}
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
MLEobj <- MARSS(dat)

library(broom)

# A data frame of the estimated parameters
tidy(MLEobj)

# Make a plot of the estimated states
library(ggplot2)
d <- tidy(MLEobj, type = "xtT")
ggplot(data = d) +
  geom_line(aes(t, estimate)) +
  geom_ribbon(aes(x = t, ymin = conf.low, ymax = conf.high), linetype = 2, alpha = 0.3) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count")

# Make a plot of the estimates for the missing values
library(ggplot2)
d <- tidy(MLEobj, type = "ytT")
ggplot(data = d) +
  geom_point(aes(t, estimate)) +
  geom_line(aes(t, estimate)) +
  geom_point(aes(t, y), color = "blue") +
  geom_ribbon(aes(x = t, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(aes(t, pred.low), linetype = 2) +
  geom_line(aes(t, pred.high), linetype = 2) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("Blue=data, Black=estimate, grey=CI, dash=prediction interval")

# Make a plot of the fitted y(t), i.e., to put a line through the points
# Intervals are for new data not the blue dots 
# (which were used to fit the model so are not new)
# Use augment() for model residuals (data - fitted values)
library(ggplot2)
d <- tidy(MLEobj, type = "fitted.ytT")
ggplot(data = d) +
  geom_line(aes(t, estimate), size=1) +
  geom_point(aes(t, y), color = "blue") +
  geom_ribbon(aes(x = t, ymin = conf.low, ymax = conf.high), alpha = 0.3) +
  geom_line(aes(t, pred.low), linetype = 2) +
  geom_line(aes(t, pred.high), linetype = 2) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("Blue=data, Black=estimate, grey=CI, dash=prediction interval") +
  geom_text(x=15, y=7, label="The intervals are for \n new data not the blue dots")
}