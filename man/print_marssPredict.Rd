\name{print.marssPredict}
\alias{print.marssPredict}

\title{ Printing function for MARSS Predict objects }
\description{
  \code{\link{MARSS}()} outputs \code{\link{marssMLE}} objects.  \code{predict(object)}, where object is \code{\link{marssMLE}} object, will return the predictions of \eqn{\mathbf{y}_t}{y(t)} or the smoothed value of \eqn{\mathbf{x}_t}{x(t)} for \code{h} steps past the end of the model data. \code{predict(object)} returns a \code{\link{marssPredict}} object which can be passed to \code{\link{print.marssPredict}} for automatic printing. 
}
\usage{
\method{print}{marssPredict}(x, ...)
}
\arguments{
  \item{x}{ A \code{\link{marssPredict}} object.  }
  \item{...}{ Other arguments for print. Not used. }
}

\value{
  A print out of the predictions as a data frame.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\examples{ 
dat <- t(harborSealWA)
dat <- dat[2:4,] #remove the year row
fit <- MARSS(dat, model=list(R="diagonal and equal"))

# 2 steps ahead forecast
predict(fit, type="ytT", h=2)

# smoothed x estimates with intervals
predict(fit, type="xtT")

}