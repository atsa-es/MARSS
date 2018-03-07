\name{augment.marssMLE}
\alias{augment.marssMLE}
\title{Return the model predicted values, residuals, and optionally confidence intervals}
\usage{
  augment.marssMLE(x, type.predict = c("observations", "states"),
                  interval = c("none", "confidence"), level = 0.95, ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type.predict}{Type of prediction: for the observations (y) or for the states (x).  Read the details below for states.  \code{tidy} would be the more common function for returning state estimates.}
  
  \item{interval}{Type of interval: none or confidence interval.  If the latter, approximate intervals from the standard errors of the fitted values is given.}
  
  \item{level}{Confidence level.}
  
  \item{...}{extra arguments (not used)}
}
\description{
Return a data.frame with the observations or states fitted values, residuals, and upper and lower confidence intervals if requested.
}
\details{

See \code{\link{residuals.marssMLE}} for a discussion of the residuals calculations for MARSS models.

\strong{type.predict observations}

This returns a familiar model predicted value of the response (y) and the difference between the model prection and the actual data y(t) is the residual.  The standard errors help visualize how well the model fits to the data.  See \code{\link{fitted.marssMLE}} for a discussion of the calculation of the fitted values for the observations (the modeled values).  The standardized residuals can be used for outlier detection.  See \code{\link{residuals.marssMLE}} and the chapter on shock detection in the user guide.

In the literature on state-space models, it is very common to use the one-step ahead predicted values of the data.  The fitted values returned by augment are NOT the one-step ahead values and the residuals are not the one-step ahead residuals (called Innovations in the state-space literature).  If you want the one-step ahead fitted values, you can use \code{fitted(x, one.step.ahead=TRUE)}.  The innovations are returned by \code{\link{MARSSkf}}.

\strong{type.predict states}

The states are estimated.  If you want the expected value of the states and an estimate of their standard errors (for confidence intervals), then \code{augment} is not what you want to use.  You want to use \code{tidy} to return the estimate of the state.

\code{augment(MLEobj, type.predict="states")} returns a model prediction of \eqn{\hat{x}(t)} given \eqn{\hat{x}(t-1)}.  The residuals returned are for w(t), the difference between \eqn{\hat{x}(t)} and the prediction of \eqn{\hat{x}(t)}.  These types of residuals are used for outlier detection or shock detection in the state process.  They are also used for model diagnostics.  See \code{\link{residuals.marssMLE}} and read the references cited.

}
