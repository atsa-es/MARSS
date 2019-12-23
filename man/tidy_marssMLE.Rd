\name{tidy.marssMLE}
\alias{tidy.marssMLE}
\title{Return estimated parameters, expected value of X(t) and Y(t) with summary information}
\description{

\code{tidy.marssMLE} returns summary information about the model parameters and estimated state and observation processes. 

For parameters, \code{tidy.marssMLE} returns their estimates and their confidence intervals. For states (X) and observations (Y), it returns the expected values (mean value) and intervals for the complete right side of the MARSS X and Y equations (including the error terms) conditioned on all the data, data up to t or data up to t-1. Do not plot the observed data on the intervals from \code{type="new observations"} as those interval prediction intervals. If you want to plot your data on intervals (e.g., for analysis of your model residuals), you need the residuals intervals. These are the intervals for model residuals (data - model fitted value). Use \code{\link{augment.marssMLE}} for the model residuals and their intervals.


The tidy function is designed to work with the \code{broom} package and you will need to load that package if you want to call \code{tidy(fit)} instead of \code{tidy.marssMLE(fit)}. 
}
\usage{
  tidy.marssMLE(x, type = c("parameters", "states", "observations", "x", "y"),
                conf.int = TRUE,
                conf.level = 0.95,
                conditioning = c("T", "t-1", "t"),
                form=attr(x[["model"]], "form")[1], ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type}{What you want estimates and intervals for. Parameters, states (x), observations (y), or new observations (y').}
  
  \item{conf.int}{ Whether to compute intervals on the estimates. }
  
  \item{conf.level}{ Confidence level. alpha=1-conf.level}

  \item{conditioning}{ If type is states or observations, how should the estimates be computed: using all data 1:T ("T"), using only data up to t-1 ("t-1"), or using only data up to t ("t"). If conditioning="t-1", then you get the one-step-ahead predictions (and their quantiles). If conditioning="T", then you get the smoothed state and model fits (and their quantiles). If type="observations", then only the values "T" and "t-1" are allowed.}

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
  
\strong{type="states"}

\code{tidy.marssMLE} returns the confidence intervals of the state at time t using the estimated model parameters as true values. These are the standard intervals that are shown for the estimated states in state-space models. For example see, Shumway and Stoffer (2000), edition 4, Figure 6.4. As such, this is probably what you are looking for if you want to put intervals on the estimated states (the x). However, these intervals do not include parameter uncertainty. If you want state residiuals (for residuals analysis), use \code{\link{residuals.marssMLE}} or \code{\link{augment.marssMLE}}.


\strong{Quantiles} The state \eqn{\mathbf{X}_t} in a MARSS model has a conditional multivariate normal distribution, that can be computed from the model parameters and data. In Holmes (2012, Eqn. 11) notation, its expected value conditioned on all the observed data (1:T) and the model parameters \eqn{\Theta} is \eqn{\tilde{\mathbf{x}}_t}{\tilde{x}_t}. In \code{\link{MARSSkf}}, this is \code{xtT[,t]}. The variance of \eqn{\mathbf{X}_t}{X_t} conditioned on the observed data and \eqn{\Theta} is \eqn{\tilde{\mathbf{V}}_t}{\tilde{V}_t} (\code{VtT[,,t]}). 

\eqn{\tilde{\mathbf{x}}_t}{\tilde{x}_t} (xtT) is an estimate of \eqn{\mathbf{x}_t}{x_t} (the true value), and the standard error of that estimate is given by \eqn{\tilde{\mathbf{V}}_t}{\tilde{V}_t} (\code{VtT[,,t]}). Let \code{se.xt} denote the sqrt of the diagonal of \code{VtT}. The equation for the \eqn{\alpha/2} confidence interval is (\code{qnorm(alpha/2)*se.xt + xtT}). \eqn{\mathbf{x}_t}{x_t} is multivariate and this interval is for one of the \eqn{x}'s in isolation. You could compute the m-dimensional confidence region for the multivariate \eqn{\mathbf{x}_t}{x_t}, also, but \code{tidy.marssMLE} returns the univariate confidence intervals. 

Analogously, we can compute the expected value and variance of \eqn{\mathbf{X}_t}{X_t} conditioned on only the data up to t-1. Those values are \code{xtt1[,t]} and \code{Vtt1[,,t]}. We can also compute the expected value and variance of \eqn{\mathbf{X}_t}{X_t} conditioned on only the data up to t. Those are \code{xtt[,t]} and \code{Vtt[,,t]}. All these statistics are available from the Kalman filter and smoother function \code{\link{MARSSkf}}.

The variances \code{VtT}, \code{Vtt} and \code{Vtt1} give information on the uncertainty of the true location of \eqn{\mathbf{x}_t} conditioned on the observed data. As more data at time t is collected (or added to the analysis), these variances will shrink since the data, especially data at time t, increases the information about the locations of \eqn{\mathbf{x}_t}. This does not affect the estimation of the model parameters, those are fixed, but rather our information about the states at time t.

If you have a DFA model (form='dfa'), you can pass in \code{rotate=TRUE} to return the rotated trends.  If you want the rotated loadings, you will need to compute those yourself:
\preformatted{
dfa <- MARSS(t(harborSealWA[,-1]), model=list(m=2), form="dfa")
Z.est <- coef(dfa, type="matrix")$Z
H.inv <- varimax(coef(dfa, type="matrix")$Z)$rotmat
Z.rot <- Z.est \%*\% H.inv
}

\strong{Intervals for the observation process}

For observation process, the expected values and intervals are shown for either new data (\code{type="new observations"}) or the observed data set (\code{type="observations"}). Details on these are below after this discussion of intervals for the observation process.

The types of intervals you want for data (Y part of the MARSS equation) depends on what you are trying to do.
\itemize{
\item{Get the model predictions}{Use \code{\link{fitted.marssMLE}} with \code{type="observations"}. This returns the fitted values (model predictions) for \eqn{Y_t} and you can condition on all the data, data up to t, or data up to t-1. Depending on what you want to do with those fitted values, you might use \code{tidy.marssMLE} or \code{\link{augment.marssMLE}}.}
\item{Get the distribution of new data that would be generated by the process}{Use \code{tidy.marssMLE} with \code{type="new observations"}.}
\item{Compare your data to model predictions}{In this case, you want the distribution of the model residuals for the data. Use \code{\link{augment.marssMLE}} with \code{type="observations"}. You want the standard errors for the data minus the fitted values (expected value of \eqn{V_t}) which is what  \code{\link{augment.marssMLE}} gives.  }
\item{Get estimates and variance of missing data in your data set}{Use \code{tidy.marssMLE} with \code{type="observations"} and \code{conditioning="T"}. The observed data will have expected value equal to the observed data and variance of 0, while the missing data will have their expected value and variance conditioned on all the observed data.}
\item{Do a leave-one-out cross-validation}{In this case, you want the distribution of the model residuals for those left-out values. Use \code{augment.marssMLE} with \code{type="observations"} (and probably \code{conditioning="T"}). You want the standard errors for the left-out data minus the fitted values (so expected value of \eqn{V_t}) which is what  \code{\link{augment.marssMLE}} gives.  }
\item{One-step-ahead predictions}{Set \code{conditioning="t-1"} for these. Both \code{\link{tidy.marssMLE}} (with \code{type="observations"} or \code{type="new observations"}) and \code{\link{augment.marssMLE}} with \code{type="observations"} give the same values in this case.  }
}

\strong{type="new observations"}

For \code{type="new observations"}, \code{tidy.marssMLE} returns the analogous information for the Y part of the MARSS equation for an I.I.D. NEW DATA SET  \eqn{y'}. The expected value and variance of \eqn{y'} is conditioned on the data you did observe \eqn{y}. It is important to note that \eqn{y'} is independent and identical (meaning i.i.d. in a statistical sense) to \eqn{y} except it has no missing values. Do not plot your observed data on these intervals. You need residuals intervals in that case. See \code{\link{augment.marssMLE}} for those. 

The expected value of a new data set \eqn{\mathbf{Y'}_t}{Y'_t} conditioned on the observed \eqn{\mathbf{Y}_t=\mathbf{y_t}}{Y_t=y_t}, is \eqn{Z_t \hat{x}_t + D_t d_t + a_t}{Z_t hatx_t + D_t d_t + a_t}, where \eqn{\hat{x}_t}{hatx_t} is the expected value of \eqn{\mathbf{X}_t}{X_t} conditioned on the data up to T, t or t-1 (depending on the value of \code{conditioning}). See \code{\link{fitted.marssMLE}}. The variance of a new data set \eqn{\mathbf{Y'}_t}{Y'_t} conditioned on the observed \eqn{\mathbf{Y}_t=\mathbf{y_t}}{Y_t=y_t}, is \eqn{R_t + Z_t \hat{V}_t Z_t^\top}{R_t + Z_t hatV_t t(Z_t)}, where \eqn{\hat{V}_t}{hatV_t} is the variance of \eqn{\mathbf{X}_t}{X_t} conditioned on the data up to T, t or t-1.

We compute the confidence interval for y', an interval that will cover the new data for alpha/2 percent of new data sets.  The equation for the \eqn{\alpha/2} confidence interval is (\code{qnorm(alpha/2)*se.yt + Ey}), where \code{Ey} is the expected value of \eqn{\mathbf{Y'}_t}{Y'_t} and \code{se.yt} is the square root of the diagonal of the variance of \eqn{\mathbf{Y'}_t}{Y'_t}.

\strong{type="observations" for missing data estimation}

This returns the expected value and variance of \eqn{\mathbf{Y}_t}{Y_t} (the data set you DID observe) conditioned on \eqn{\mathbf{Y}_t=y_t}{Y_t=y_t}. This what you want if you want to estimate the values of missing data in your data set. The expected value of \eqn{\mathbf{Y}_t|\mathbf{Y}_t=\mathbf{y}_t}{Y_t|Y_t=y_t} is in \code{ytT} in \code{\link{MARSShatyt}} output and the variance is \code{OtT-tcrossprod(ytT)} from the \code{\link{MARSShatyt}} output. 

The expected value of \eqn{\mathbf{Y}_t}{Y_t} conditioned on \eqn{\mathbf{Y}_t=\mathbf{y}_t}{Y_t=y_t} is simply \eqn{\mathbf{y}_t}{y_t}. So is not useful if you have no missing values. Its use is for the case where you have missing values scattered throughout your data set. The standard error and confidence interval reported by \code{tide.marssMLE} for the missing values takes into account all the information in the data, specifically the correlation with other data at time t if R is not diagonal. Do not use \code{\link{augment.marssMLE}} as that returns the standard errors for the model residuals for missing data, which will be quite different if the R matrix is not diagonal (though identical if R is diagonal). Also do not use "prediction" intervals for interpolating missing data (\code{type="new observations"}), as those are for entirely new data sets and thus will ignore relevant information if \eqn{\mathbf{y}_t}{y_t} is multivariate, not all \eqn{\mathbf{y}_t}{y_t} are missing, and the R matrix is not diagonal.

\strong{One step ahead predictions} Use \code{conditioning="t-1"} for the one-step-ahead predictions of \eqn{\mathbf{y}_t}{y_t}, using only data up to t-1. Because the one-step-ahead predictions do not include \eqn{y_t} for calculation, you can use your data to evaluate the intervals. Thus in this case (only) it would make sense to plot your data on top of the prediction intervals. But in general use \code{\link{augment.marssMLE}} when you want statistics (such as intervals) associated with your observed data.

\strong{Output that \code{augment.marssMLE} provides and \code{tidy.marssMLE} does not}


\strong{fitted Y(t)} The fitted values of \eqn{\mathbf{Y}_t}{Y_t} is the expected value of the right side of the \eqn{\mathbf{Y}_t}{Y_t} equation WITHOUT the \eqn{v_t}. Use \code{\link{fitted.marssMLE}} and \code{\link{augment.marssMLE}} to get the fitted values, their residuals, and their residual variances (for outlier analysis).

\strong{Residual intervals} For residuals analysis (outlier detection), you want to use the "residual" intervals provided with \code{\link{augment.marssMLE}}. These are are the quantiles for the model residuals for data that were used to fit the model (and estimate xtT).  Read more about residuals analysis for MARSS models at \code{\link{residuals.marssMLE}}. The residuals intervals are also what you want to evaluate data that you leave out, for example for a leave-one-out-cross-validation.

\strong{Leave-one-out analysis} The one-step-ahead prediction intervals are not what you want if you left some of your data out and want to evaluate the model residuals for the left-out-data. The residuals intervals are also what you want for this case (with the left-out-data set to NA).

\strong{Prediction intervals for new data} The one-step-ahead prediction intervals are also not what you want if you want to compute intervals for data during 1 to TT that you didn't collect (so not in the future but contemporaneous with your data).  To compute prediction intervals for new data within 1:TT (so not forecasts past TT), you need to decide what you mean by 'new' and how these new data at time t relate to existing data at time t. Are the new data are i.i.d relative to the old data? So if you made an R matrix, it would be diagonal with the same value on the diagonal. Or are they missing data which might have a different variance and be correlated with the old data. In this case the R matrix might not be diagonal with one value on the diagonal.

\strong{Parameter uncertainty} Currently the intervals calculations for the states and observations use the point estimates of the model parameters and thus solve the intervals for the 'known' parameters case.

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
  facet_grid(~.rownames) +
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
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("Blue=data, Black=estimate of missing data point")
  
}
