\name{tidy.marssMLE}
\alias{tidy.marssMLE}
\title{Return estimated parameters, expected value of X(t) and Y(t) with summary information}
\description{

For parameters, \code{tidy.marssMLE} returns their estimates and their confidence intervals. For states (X) and observations (Y), it returns the expected values (mean value) and intervals for the complete right side of the MARSS X and Y equations (including the error terms) conditioned on all the data, data up to t or data up to t-1. 

If you want the fitted values of \eqn{\mathbf{Y}_t}{Y_t}, use \code{\link{fitted.marssMLE}} or \code{\link{augment.marssMLE}}. \code{\link{augment.marssMLE}} is  also what you want for residuals analysis (for outlier analysis). See details.

The tidy function is designed to work with the \code{broom} package and you will need to load that package if you want to call \code{tidy(fit)} instead of \code{tidy.marssMLE(fit)}. 
}
\usage{
  tidy.marssMLE(x, type = c("parameters", "states", "observations", "x", "y"),
                conf.int = TRUE,
                conf.level = 0.95,
                smoothing = c("T", "t-1", "t"),
                form=attr(x[["model"]], "form")[1], ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type}{What you want estimates and intervals for. Parameters, states (X) or observations (Y).}
  
  \item{conf.int}{ Whether to compute intervals on the estimates. }
  
  \item{conf.level}{ Confidence level. alpha=1-conf.level}

  \item{smoothing}{ If type is states or observations, how should the estimates be computed: using all data 1:T ("T"), using only data up to t-1 ("t-1"), or using only data up to t ("t"). If smoothing="t-1", then you get the one-step-ahead predictions (and their quantiles). If smoothing="T", then you get the smoothed state and model fits (and their quantiles). }

  \item{form}{If you want the tidy function to use a different form than that specified in attr(x$model, "form").  Useful if you have a DFA model that you manually set up, which does not have the form attribute set. Normally just ignore and let the function use the "form" set in the attributes.}
  
  \item{...}{Optional arguments.  If \code{conf.int=TRUE}, then arguments to specify how CIs are computed can be passed in. See details and \code{\link{MARSSparamCIs}}. If form="dfa", \code{rotate=TRUE} can be passed in to rotate the trends (only trends not Z matrix).}

}

\value{
A data.frame with estimates, sample standard errors, and confidence (or prediction) intervals.
}

\details{

Below, X and Y refers to the random variable and x and y refer to the actual value (which is drawn from the distribution of the random variable).

\strong{type="parameters"}  

If \code{type="parameters"}, this returns a data.frame with the estimated parameters  of a MARSS model with, optionally, standard errors and confidence intervals. This assembles information available via the \code{\link{print.marssMLE}} and \code{\link{coef.marssMLE}} functions into a data.frame that summarizes the estimates.  If conf.int=TRUE, \code{\link{MARSSparamCIs}} will be run to add confidence intervals to the fitted model object if these are not already added to the model object.  The default CIs are calculated using a analytically computed Hessian matrix.  This can be changed by passing in optional arguments for \code{\link{MARSSparamCIs}}. 
  
\strong{type="states"}

\code{tidy.marssMLE} returns the confidence intervals of the state at time t using the estimated model parameters as true values. These are the standard intervals that are shown for the estimated states in state-space models. For example see, Shumway and Stoffer (2000), edition 4, Figure 6.4. As such, this is probably what you are looking for if you want to put intervals on the estimated states (the x). However, these intervals do not include parameter uncertainty. If you want state residiuals (for residuals analysis), use \code{\link{residuals.marssMLE}} or \code{\link{augment.marssMLE}}.


\strong{Quantiles} The state \eqn{\mathbf{X}_t} in a MARSS model has a conditional multivariate normal distribution, that can be computed from the model parameters and data. In Holmes (2012, Eqn. 11) notation, its expected value conditioned on all the observed data (1:T) and the model parameters \eqn{\Theta} is \eqn{\tilde{\mathbf{x}}_t}{\tilde{x}_t}. In \code{\link{MARSSkf}}, this is \code{xtT[,t]}. The variance of \eqn{\mathbf{X}_t}{X_t} conditioned on the observed data and \eqn{\Theta} is \eqn{\tilde{\mathbf{V}}_t}{\tilde{V}_t} (\code{VtT[,,t]}). 

\eqn{\tilde{\mathbf{x}}_t}{\tilde{x}_t} (xtT) is an estimate of \eqn{\mathbf{x}_t}{x_t} (the true value), and the standard error of that estimate is given by \eqn{\tilde{\mathbf{V}}_t}{\tilde{V}_t} (\code{VtT[,,t]}). Let \code{se.xt} denote the sqrt of the diagonal of \code{VtT}. The equation for the \eqn{\alpha/2} confidence interval is (\code{qnorm(alpha/2)*se.xt + xtT}). \eqn{\mathbf{x}_t}{x_t} is multivariate and this interval is for one of the \eqn{x}'s in isolation. You could compute the m-dimensional confidence region for the multivariate \eqn{\mathbf{x}_t}{x_t}, also, but \code{tidy.marssMLE} returns the univariate confidence intervals. 

Analogously, we can compute the expected value and variance of \eqn{\mathbf{X}_t}{X_t} conditioned on only the data up to t-1. Those values are \code{xtt1[,t]} and \code{Vtt1[,,t]}. We can also compute the expected value and variance of \eqn{\mathbf{X}_t}{X_t} conditioned on only the data up to t. Those are \code{xtt[,t]} and \code{Vtt[,,t]}. All these statistics are available from the Kalman filter and smoother function \code{\link{MARSSkf}}.

If you have a DFA model (form='dfa'), you can pass in \code{rotate=TRUE} to return the rotated trends.  If you want the rotated loadings, you will need to compute those yourself:
\preformatted{
dfa <- MARSS(t(harborSealWA[,-1]), model=list(m=2), form="dfa")
Z.est <- coef(dfa, type="matrix")$Z
H.inv <- varimax(coef(dfa, type="matrix")$Z)$rotmat
Z.rot <- Z.est \%*\% H.inv
}


\strong{type="observations"}

\code{tidy.marssMLE} returns the expected value and variance of \eqn{\mathbf{Y}_t}{Y_t}, conditioned on all the data (conditioning on data up to t or t-1 is not provided). Note that the expected value of \eqn{\mathbf{Y}_t}{Y_t} conditioned on \eqn{\mathbf{Y}_t=\mathbf{y_t}}{Y_t=y_t}, meaning if there is no missing data, the expected value of the data is the data. However, if you have missing data, the expected value and variance of \eqn{\mathbf{Y}_t}{Y_t} is what you want.  These values are returned by \code{\link{MARSShatyt}} in \code{ytT} and \code{OtT-tcrossprod(ytT)}. The equation for the \eqn{\alpha/2} confidence interval is (\code{qnorm(alpha/2)*se.yt + Ey}). Note that \code{ytT} is not the same as the fitted \eqn{\mathbf{Y}_t}{Y_t}, and if there are no missing values, the variance of \eqn{\mathbf{Y}_t}{Y_t} conditioned on the data will be 0. See below.  

\strong{The estimate and distribution for missing data} If your \eqn{\mathbf{y}_t}{y_t} is missing data and you want an estimate of those values, use \code{tidy.marssMLE}. The standard error and confidence interval for the missing value takes into account all the information in the data. Do not use \code{\link{augment.marssMLE}} as that returns the standard errors for the model residuals for missing data, which will be quite different if the R matrix is not diagonal (though identical if R is diagonal). Also do not use "prediction" intervals for interpolating missing data, as those are for entire new data sets and thus will ignore relevant information if \eqn{\mathbf{y}_t}{y_t} is multivariate, not all \eqn{\mathbf{y}_t}{y_t} are missing, and the R matrix is not diagonal.

\strong{Output that \code{augment.marssMLE} provides and \code{tidy.marssMLE} does not}

\strong{fitted Y(t)} The fitted values of \eqn{\mathbf{Y}_t}{Y_t} is the expected value of the right side of the \eqn{\mathbf{Y}_t}{Y_t} equation WITHOUT the \eqn{v_t}. Use \code{\link{fitted.marssMLE}} and \code{\link{augment.marssMLE}} to get the fitted values, their residuals, and their residual variances (for outlier analysis).

\strong{Residual intervals} For residuals analysis (outlier detection), you want to use the "residual" intervals provided with \code{\link{augment.marssMLE}}. These are are the quantiles for the model residuals for data that were used to fit the model (and estimate xtT).  Read more about residuals analysis for MARSS models at \code{\link{residuals.marssMLE}}. The residuals intervals are also what you want to evaluate data that you leave out, for example for a leave-one-out-cross-validation.

\strong{Prediction intervals} The prediction interval is for new independent data which observe the same \eqn{\mathbf{x}_t}{x_t} as the data which were used to fit the model. Use \code{\link{augment.marssMLE}} for prediction intervals.

Currently the intervals calculations for the states and observations use the point estimates of the model parameters and thus solve the intervals for the 'known' parameters case.
}
\references{ 
R. H. Shumway and D. S. Stoffer (2000).  Time series analysis and its applications. Edition 4. Springer-Verlag, New York.

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]
}
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11,12),]
  # model allowing correlation between x states
  MLEobj <- MARSS(dat, model=list(Q="unconstrained"))

library(broom)
library(ggplot2)

# A data frame of the estimated parameters
 tidy(MLEobj)
 
# Make a plot of the estimated states
library(ggplot2)
d <- tidy(MLEobj, type="states")
ggplot(data = d) + 
  geom_line(aes(t, estimate)) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.3) +
  facet_grid(~term) +
  xlab("Time Step") + ylab("Count")
  
# Make a plot of the estimates for the missing values
# Use augment() to plot the fitted y(t)
# i.e., to put a line through the points
library(ggplot2)
d <- tidy(MLEobj, type="observations")
ggplot(data = d) + 
  geom_point(aes(t, estimate)) +
  geom_line(aes(t, estimate)) +
  geom_point(aes(t, y), color="blue") +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.3) +
  facet_grid(~term) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("Blue=data, Black=estimate of missing data point")
  
}
