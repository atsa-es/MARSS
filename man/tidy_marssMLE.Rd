\name{tidy.marssMLE}
\alias{tidy.marssMLE}
\title{Return estimated parameters and states with summary information}
\usage{
  tidy.marssMLE(x, type = c("parameters", "states"), conf.int = FALSE,  ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type.predict}{Type of prediction: for the observations (y) or for the states (x).  Read the details below for states.  \code{tidy} would be the more common function for returning state estimates.}
  
  \item{conf.int}{ Whether to include a confidence interval.}
  
  \item{...}{optional arguments to be passed to \code{\link{MARSSparamCIs}} if conf.int is TRUE.  See \code{\link{MARSSparamCIs}} for the default values and details.}:
  \item{alpha}{Confidence level}
  
  \item{method}{Method to use to calculate the confidence interval.  Analytical and numerical options are available.}
  
  \item{hessian.fun}{Function to use to compute the Hessian if method="hessian".}
  
  \item{nboot}{number of bootstraps for the parametric or innovations methods.}
}
\description{
  This returns a data.frame with the estimated parameters (or states) of a MARSS model with standard errors and optionally confidence intervals. This assembles information available via the print and coef functions into a data.frame that summarizes the estimates and can be used for further processing.
  
  The tidy function is compatible with the \code{broom} package.
}
\details{
}