\name{glance}
\alias{glance}
\alias{glance.marssMLE}
\keyword{coremethods}
\title{Return brief summary information on a MARSS fit}
\usage{
\method{glance}{marssMLE}(x, ...)
}
\arguments{
  \item{x}{A \code{\link{marssMLE}} object}

  \item{...}{Not used.}
}
\description{
  This returns a data frame with brief summary information.
  \describe{
    \item{coef.det}{ The coefficient of determination. This is the squared correlation between the fitted values and the original data points. This is simply a metric for the difference between the data points and the fitted values and should not be used for formal model comparison.  }
    \item{sigma}{ The sample variance (unbiased) of the data residuals (fitted minus data). This is another simple metric of the difference between the data and fitted values.  This is different than the sigma returned by an \code{\link[stats]{arima}()} call for example.  That sigma would be akin to \code{sqrt(Q)} in the MARSS parameters; 'akin' because MARSS models are multivariate and the sigma returned by \code{arima()} is for a univariate model. }
    \item{df}{ The number of estimated parameters.  Denoted \code{num.params} in a \code{\link{marssMLE}} object.}
    \item{logLik}{ The log-likelihood.}
    \item{AIC}{ Akaike information criterion.}
    \item{AICc}{ Akaike information criterion corrected for small sample size.}
    \item{AICbb}{ Non-parametric bootstrap Akaike information criterion if in the \code{\link{marssMLE}} object.}
    \item{AICbp}{ Parametric bootstrap Akaike information criterion if in the \code{\link{marssMLE}} object.}
    \item{convergence}{ 0 if converged according to the convergence criteria set.  Note the default convergence criteria are high in order to speed up fitting. A number other than 0 means the model did not meet the convergence criteria.}
    \item{errors}{ 0 if no errors. 1 if some type of error or warning returned.}
  }
}
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
fit <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))

glance(fit)
}