\name{MARSShessian.numerical}
\alias{MARSShessian.numerical}
\title{ MARSS Hessian Matrix via Numerical Approximation}
\description{
Calculates the Hessian of the log-likelihood function at the MLEs using either the  \code{\link{fdHess}} function in the  \code{\link{nmle}} package or the  \code{\link{optim}} function.  This is a utility function in the \code{\link{MARSS-package}} and is not exported.  Use \code{\link{MARSShessian}} to access.
}
\usage{
MARSShessian.numerical(MLEobj, fun="fdHess")
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. 
  This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{\link{MARSSkem}}. }
  \item{fun}{ The function to use for computing the Hessian. Options are 'fdHess' or 'optim'.}
 }
\details{
   Method 'fdHess' uses \code{\link{fdHess}} from package \code{\link{nlme}} to numerically estimate the Hessian matrix
   (the matrix of partial 2nd derivatives of the parameter estimates).     Method 'optim' uses \code{\link{optim}} with \code{hessian=TRUE} and \code{list(maxit=0)} to ensure that the Hessian is computed at the values in the \code{par} element of the MLE object.
}
\value{
  The numerically estimated Hessian of the log-likelihood function at the MLEs.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSSharveyobsFI}} \code{\link{MARSShessian}} \code{\link{MARSSparamCIs}} }
\examples{
  dat = t(harborSeal)
  dat = dat[c(2,11),]
  MLEobj = MARSS(dat)
  MARSS:::MARSShessian.numerical(MLEobj)
}

