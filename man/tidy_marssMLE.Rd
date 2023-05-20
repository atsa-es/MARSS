\name{tidy.marssMLE}
\alias{tidy}
\alias{tidy.marssMLE}
\title{Return estimated parameters with summary information}
\description{

\code{tidy.marssMLE} is the method for the tidy generic. It returns  the parameter estimates and their confidence intervals. 
}
\usage{
\method{tidy}{marssMLE}(x, conf.int = TRUE, conf.level = 0.95, ...)
}
\arguments{
  \item{x}{a \code{\link{marssMLE}} object}
  
  \item{conf.int}{ Whether to compute confidence and prediction intervals on the estimates. }
  
  \item{conf.level}{ Confidence level. \code{alpha=1-conf.level}}

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
tidy(fit)

}