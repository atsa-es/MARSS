\name{stdInnov}
\alias{stdInnov}
\keyword{internal}

\title{ Standardized Innovations }
\description{
  Standardizes Kalman filter innovations. This is a helper function called by \code{\link{MARSSinnovationsboot}()} in the \code{\link{MARSS-package}}. Not exported.
}
\usage{
stdInnov(SIGMA, INNOV)
}
\arguments{
  \item{SIGMA}{ n x n x T array of Kalman filter innovations variances.  This is output from \code{\link{MARSSkf}}. }
  \item{INNOV}{ n x T matrix of Kalman filter innovations.  This is output from  \code{\link{MARSSkf}()}. }
}
\details{
  n = number of observation (y) time series.   T = number of time steps in the time series.
}
\value{
  n x T matrix of standardized innovations.
}
\references{ 
Stoffer, D. S., and K. D. Wall. 1991. Bootstrapping state-space models: Gaussian maximum likelihood estimation and the Kalman filter. Journal of the American Statistical Association 86:1024-1033.
 }

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\seealso{  \code{\link{MARSSboot}}, \code{\link{MARSSkf}}, \code{\link{MARSSinnovationsboot}}
}
\examples{
  \dontrun{
  std.innovations <- stdInnov(kfList$Sigma, kfList$Innov)
  }
}

