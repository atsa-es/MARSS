\name{zscore}
\alias{zscore}
\keyword{user-helper}
\title{ z-score a vector or matrix }
\description{
  Removes the mean and standardizes the variance to 1.
}
\usage{
zscore(x, mean.only = FALSE)
}
\arguments{
  \item{x}{ n x T matrix of numbers }
  \item{mean.only}{ If TRUE, only remove the mean. }
}
\details{
  n = number of observation (y) time series.   T = number of time steps in the time series.
  
  The z-scored values (z) of a matrix of y values are
\eqn{z_i = \Sigma^{-1}(y_i-\bar{y})}
  where \deqn{\Sigma} is a diagonal matrix with the standard deviations of each time series (row) along the diagonal, and \deqn{\bar{y}} is a vector of the means. 
}
\value{
  n x T matrix of z-scored values.
}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\examples{
zscore(1:10)
x <- zscore(matrix(c(NA, rnorm(28), NA), 3, 10))
# mean is 0 and variance is 1
apply(x, 1, mean, na.rm = TRUE)
apply(x, 1, var, na.rm = TRUE)
}

