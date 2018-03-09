\name{tidy.marssMLE}
\alias{tidy.marssMLE}
\title{Return estimated parameters and states with summary information}
\usage{
  tidy.marssMLE(x, type = c("parameters", "states"), conf.int = TRUE,  ...)
}
\arguments{
  \item{x}{a marssMLE object}
  
  \item{type}{Estimates for the parameters or for the states.}
  
  \item{conf.int}{ Whether to include a confidence interval.}
  
  \item{...}{optional arguments to be passed to \code{\link{MARSSparamCIs}} if conf.int is TRUE.  See \code{\link{MARSSparamCIs}} for the default values and details.}:

  \item{alpha}{Confidence level}
  
  \item{method}{Method to use to calculate the confidence interval.  Analytical and numerical options are available.}
  
  \item{hessian.fun}{Function to use to compute the Hessian if method="hessian".}
  
  \item{nboot}{number of bootstraps for the parametric or innovations methods.}
}
\description{
  This returns a data.frame with the estimated parameters (or states) of a MARSS model with optionally standard errors and confidence intervals. This assembles information available via the print and coef functions into a data.frame that summarizes the estimates.
  
  For parameters, \code{\link{MARSSparamsCIs}} is run to add confidence intervals to the fitted model if these are not already added.  For states, the approximate CIs using the standard deviation of the states is used to compute the confidence intervals (\code{qnorm(alpha/2)*se.fit + fitted}).
  
  The tidy function is compatible with the \code{broom} package.
}
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11,12),]
  MLEobj <- MARSS(dat, model=list(Z=factor(c("WA","OR","OR"))))

library(broom)

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