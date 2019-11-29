\name{tidy.marssMLE}
\alias{tidy.marssMLE}
\title{Return estimated parameters and states with summary information}
\usage{
  tidy.marssMLE(x, type = c("parameters", "states"),
                interval = TRUE, alpha = 0.05,
                form=attr(x[["model"]], "form")[1], ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type}{Estimates for the parameters or for the states.}
  
  \item{interval}{ Whether to include a confidence interval if type="parameters" and quantiles if type is "states", "observations", or any of the x or y variables. See details for discussion of intervals for non-parameters.}

  \item{alpha}{ Confidence level if confidence interval is returned otherwise the quantile level for quantiles.}

  \item{form}{If you want the augment function to use a different form than that specified in attr(x$model, "form").  Useful if you have a DFA model that you manually set up, which does not have the form attribute set.}
  
  \item{...}{Optional arguments.  If \code{interval=TRUE}, then arguments to specify how CIs are computed can be passed in. See details and \code{\link{MARSSparamCIs}}. If form="dfa", \code{rotate=TRUE} can be passed in to rotate the trends (only trends not Z matrix).}

}
\description{

The tidy function is designed to work with the \code{broom} package and you will need to load that package if you want to call \code{tidy(fit)} instead of \code{tidy.marssMLE(fit)}.

\strong{type=="parameters"}

  If type=="parameters", this returns a data.frame with the estimated parameters  of a MARSS model with, optionally, standard errors and confidence intervals. This assembles information available via the print and coef functions into a data.frame that summarizes the estimates.
  
  If interval=TRUE, \code{\link{MARSSparamCIs}} will be run to add confidence intervals to the fitted model object if these are not already added to the model object.  The default CIs are calculated using a analytically computed Hessian matrix.  This can be changed by passing in optional arguments for \code{\link{MARSSparamCIs}}. 
  
\strong{type=="states"}

Currently \code{tidy.marssMLE} does not return the confidence intervals for the estimated states. Instead only the quantiles of the states at time (t) using the estimated model parameters is returned. These are the standard intervals that are shown for the estimated states in state-space models. For example see, Shumway and Stoffer (2000), edition 4, Figure 6.4. As such, this is probably what you are looking for if you want to put intervals on the estimated states (the X). However, see the comments below on confidence intervals for the expected value of the states.

\strong{Quantiles} The state X(t) in a MARSS model is a random variable. Using \code{MARSSkf} notation, the estimate of its expected value conditioned on all the observed data is xtT(t). It's an estimate because the Theta we use is an estimate, i.e. \hat{Theta}. The estimate of the variance of X(t) conditioned on all the observed data and \hat{Theta}, i.e. E[tcrossprod(X(t)-xtT(t))], is VtT(t). Analogously, we can compute the expected value and variance of X(t) conditioned on only the data up to t-1. Those values are xtt1(t) and Vtt1(t). We can also compute the expected value and variance of X(t) conditioned on only the data up to t. Those are xtt(t) and Vtt(t). All these statistics are available from the Kalman filter and smoother function \code{MARSSkf}.

For states, \code{tidy.marssMLE} returns the empirical, i.e. estimated, quantiles of x(t) based on one of the estimated conditional variance-covariance matrices and corresponding expected value of X(t). Let's denote the former as Vx (which is either VtT, Vtt1 or Vtt) and the latter as Ex (which is either xtT, xtt1, or xtt). Let sd.xt denote the sqrt of the diagonal of Vx(t).  The equation for the alpha/2 quantile is (\code{qnorm(alpha/2)*sd.xt + Ex}). To emphasize that these are intervals for X(t) not E[X(t)], the column headings are .quant.low and .quant.up and the column on its standard deviation is labelled .std.dev not .std.error. A standard error is a measure of the variability of some estimated value from a sample and this variability shrinks a the sample size increases. The standard deviation does not shrink with more data; it is simply a function of the model parameters.
  
If you have a DFA model (form='dfa'), you can pass in \code{rotate=TRUE} to return the rotated trends.  If you want the rotated loadings, you will need to compute those yourself:
\preformatted{
dfa <- MARSS(t(harborSealWA[,-1]), model=list(m=2), form="dfa")
Z.est <- coef(dfa, type="matrix")$Z
H.inv <- varimax(coef(dfa, type="matrix")$Z)$rotmat
Z.rot <- Z.est \%*\% H.inv
}

\strong{Confidence intervals for \hat{E[X(t)]}} Although \code{tidy.marssMLE} does not return these currently, computation of them is conceptually straight forward. The standard error of \hat{E[X(t)]}, whether xtT, xtt1 or xtt, is a measure of the how far \hat{E[X(t)]} computed with \hat{Theta} is from the true E(X(t)) with the true Theta. The standard error of \hat{E[X(t)]} uses the information about our uncertainty in \hat{Theta}. With more data, i.e. a longer time series or multiple y(t) at each t, the uncertainty in our estimated Theta shrinks and the standard error of \hat{E[X(t)]} shrinks. This is the standard notion of the standard error of a sample statistic.

So how do you get the standard error of \hat{E[X(t)]}? 
  
}
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11,12),]
  MLEobj <- MARSS(dat, model=list(Z=factor(c("WA","OR","OR"))))

library(broom)
library(ggplot2)

# A data frame of the estimated parameters
 tidy(MLEobj)
 
# Make a plot of the estimated states
# Don't use augment.  States are not data.
d <- tidy(MLEobj, type="states")
ggplot(data = d) + 
  geom_line(aes(t, estimate)) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1) +
  facet_grid(~term) +
  xlab("Time Step") + ylab("Count")
}