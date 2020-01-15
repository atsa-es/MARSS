\name{MARSSinnovationsboot}
\alias{MARSSinnovationsboot}
\title{ Bootstrapped Data using Stoffer and Wall's Algorithm }
\description{
  Creates bootstrap data via sampling from the standardized innovations matrix.   This is an internal function in the \code{\link{MARSS-package}} and is not exported. Users should access this with \code{\link{MARSSboot}}.
}
\usage{
MARSSinnovationsboot(MLEobj, nboot = 1000, minIndx = 3)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{\link{MARSSkem}()} or \code{\link{MARSS}()}. This algorithm may not be used if there are missing datapoints in the data. }
  \item{nboot}{ Number of bootstraps to perform. }
  \item{minIndx}{ Number of innovations to skip. Stoffer & Wall suggest not sampling from innovations 1-3. }
}
\details{
  Stoffer and Wall (1991) present an algorithm for generating CIs via a non-parametric bootstrap for state-space models.  The basic idea is that the Kalman filter can be used to generate estimates of the residuals of the model fit.  These residuals are then standardized and resampled and used to generate bootstrapped data using the MARSS model and its maximum-likelihood parameter estimates.  One of the limitations of the Stoffer and Wall algorithm is that it cannot be used when there are missing data, unless all data at time t are missing.
}
\value{
  A list containing the following components:
  \item{boot.states}{ Array (dim is m x tSteps x nboot) of simulated state processes. }
  \item{boot.data}{ Array (dim is n x tSteps x nboot) of simulated data. }
  \item{marss}{ \code{\link{marssMODEL}} object element of the \code{\link{marssMLE}} object (\code{marssMLE$marss}) in "marss" form. }
  \item{nboot}{ Number of bootstraps performed. }
  m is the number state processes (x in the MARSS model) and n is the number of observation time series (y in the MARSS model).
}
\references{ 
Stoffer, D. S., and K. D. Wall. 1991. Bootstrapping state-space models: Gaussian maximum likelihood estimation and the Kalman filter. Journal of the American Statistical Association 86:1024-1033.
 }
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ \code{\link{stdInnov}}, \code{\link{MARSSparamCIs}}, \code{\link{MARSSboot}} }
\examples{
dat <- t(kestrel)
dat <- dat[2:3, ]
MLEobj <- MARSS(dat, model = list(U = "equal", Q = diag(.01, 2)))
boot.obj <- MARSSinnovationsboot(MLEobj)
}
