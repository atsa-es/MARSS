\name{augment.marssMLE}
\alias{augment.marssMLE}
\alias{augment_dfa}
\alias{augment_marss}
\alias{augmentmarxss}
\title{Return the model predicted values, residuals, and optionally confidence intervals}
\usage{
  augment.marssMLE(x, type.predict = c("observations", "states"),
                              interval = c("none", "confidence"), 
                              conf.level = 0.95, form=attr(x[["model"]], "form"), ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type.predict}{Type of prediction: for the observations (y) or for the states (x).  Read the details below for states.  \code{tidy} would be the more common function for returning state estimates.}
  
  \item{interval}{Type of interval: none or confidence interval.  If the latter, approximate intervals from the standard errors of the fitted values is given.}
  
  \item{conf.level}{Confidence level.}
  
  \item{form}{If you want the augment function to use a different augment function than augment.form.  This might be useful if you manually specified a DFA model and want to use augment.dfa for rotating.}

  \item{...}{extra arguments (not used)}
}
\description{
Return a data.frame with the observations or states fitted values, residuals, and upper and lower confidence intervals if requested.
}
\details{

See \code{\link{residuals.marssMLE}} for a discussion of the residuals calculations for MARSS models.
The reported CIs are the approximate CIs computed using the standard deviations: \code{qnorm(alpha/2)*se.fit + fitted}.

\strong{type.predict observations}

This returns a familiar model predicted value of the response (y) and the difference between the model prediction and the actual data y(t) is the residual.  If there is no data point, then the residual is NA.  The standard errors help visualize how well the model fits to the data.  See \code{\link{fitted.marssMLE}} for a discussion of the calculation of the fitted values for the observations (the modeled values).  The standardized residuals can be used for outlier detection.  See \code{\link{residuals.marssMLE}} and the chapter on shock detection in the user guide.

In the literature on state-space models, it is very common to use the one-step ahead predicted values of the data.  The fitted values returned by augment are NOT the one-step ahead values and the residuals are not the one-step ahead residuals (called Innovations in the state-space literature).  If you want the one-step ahead fitted values, you can use \code{fitted(x, one.step.ahead=TRUE)}.  The innovations are returned by \code{\link{MARSSkf}}.

\strong{type.predict states}

The states are estimated.  If you want the expected value of the states and an estimate of their standard errors (for confidence intervals), then \code{augment} is not what you want to use.  You want to use \code{tidy} to return the estimate of the state.

\code{augment(MLEobj, type.predict="states")} returns a model prediction of \eqn{\hat{x}(t)} given \eqn{\hat{x}(t-1)}.  The residuals returned are for w(t), the difference between \eqn{\hat{x}(t)} and the prediction of \eqn{\hat{x}(t)}.  These types of residuals are used for outlier detection or shock detection in the state process.  They are also used for model diagnostics.  See \code{\link{residuals.marssMLE}} and read the references cited.

}

\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11,12),]
  MLEobj <- MARSS(dat, model=list(Z=factor(c("WA","OR","OR"))))

library(broom)
library(ggplot2)
theme_set(theme_bw())

# Make a plot of the observations and model fits
d <- augment(MLEobj, interval="confidence")
ggplot(data = d) + 
  geom_line(aes(t, .fitted)) +
  geom_point(aes(t, y)) +
  geom_ribbon(aes(x=t, ymin=.conf.low, ymax=.conf.up), linetype=2, alpha=0.1) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count")

# Make a plot of the estimated states
# Don't use augment.  States are not data.
d <- tidy(MLEobj, type="states")
ggplot(data = d) + 
  geom_line(aes(t, estimate)) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1) +
  facet_grid(~term) +
  xlab("Time Step") + ylab("Count")

}