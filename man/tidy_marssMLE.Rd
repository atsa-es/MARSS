\name{tidy.marssMLE}
\alias{tidy.marssMLE}
\title{Return estimated parameters with summary information}
\description{

\code{tidy.marssMLE} returns  the parameter estimates and their confidence intervals. 

The tidy function is designed to work with the \code{broom} package and you will need to load that package if you want to call \code{tidy(fit)} instead of \code{tidy.marssMLE(fit)}. 
}
\usage{
tidy.marssMLE(x, conf.int = TRUE, conf.level = 0.95, ...)
}
\arguments{
  \item{x}{a \code{\link{marssMLE}} object}
  
  \item{conf.int}{ Whether to compute confidence and prediction intervals on the estimates. }
  
  \item{conf.level}{ Confidence level. alpha=1-conf.level}

  \item{...}{Optional arguments.  If \code{conf.int=TRUE}, then arguments to specify how CIs are computed can be passed in. See details and \code{\link{MARSSparamCIs}}. }

}

\value{
A data frame with estimates, sample standard errors, and confidence intervals.
}

\details{

\code{tidy.marssMLE()} assembles information available via the \code{\link[=print.marssMLE]{print}()} and \code{\link[=coef.marssMLE]{coef}()} functions into a data frame that summarizes the estimates.  If \code{conf.int=TRUE}, \code{\link{MARSSparamCIs}()} will be run to add confidence intervals to the model object if these are not already added.  The default CIs are calculated using a analytically computed Hessian matrix.  This can be changed by passing in optional arguments for \code{\link{MARSSparamCIs}()}. 
  
}

\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
fit <- MARSS(dat)

# A data frame of the estimated parameters
tidy.marssMLE(fit)

# If broom is used

library(broom)
tidy(fit)

}