\name{MARSShessian}
\alias{MARSShessian}
\alias{MARSSFisherI}
\title{ Parameter Variance-Covariance Matrix from the Hessian Matrix}
\description{
Calculates an approximate parameter variance-covariance matrix for the parameters using an inverse of the Hessian of the log-likelihood function at the MLEs.   It appends \code{$Hessian}, \code{$parMean}, \code{$parSigma} to the marssMLE object. 
}
\usage{
MARSShessian( MLEobj, method=c("Harvey1989", "fdHess", "optim") )
MARSSFisherI( MLEobj, method=c("Harvey1989", "fdHess", "optim") )
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. 
  This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{\link{MARSSkem}}. }
  \item{method}{ The method to use for computing the Hessian. Options are 'Harvey1989' to use the Harvey (1989) recursion, which is an analytical solution, 'fdHess' or 'optim' which are two numerical methods.  Although 'optim' can be passed to the function, 'fdHess' is used for all numerical estimates used in the package.}
 }
\details{
   Method 'fdHess' uses \code{\link{fdHess}} from package \code{\link{nlme}} to numerically estimate the Hessian matrix
   (the matrix of partial 2nd derivatives of the parameter estimates).     Method 'optim' uses \code{\link{optim}} with \code{hessian=TRUE} and \code{list(maxit=0)} to ensure that the Hessian is computed at the values in the \code{par} element of the MLE object.
   Method 'Harvey1989' (the default) uses the recursion in Harvey (1989) to compute the observed Fisher Information of a MARSS model analytically.
   
   Hessian CIs are based on 
   the asymptotic normality of ML estimates under a large-sample approximation. 
   
   MARSSFisherI() is a (non-exported) utility function called by MARSShessian().  It calls either \code{\link{MARSShessian.numerical}} or \code{\link{MARSSharveyobsFI}} to return the Hessian based on the value of \code{method}.
   
}
\value{
  \code{MARSShessian()} attaches 
   \code{Hessian}, \code{parMean} and \code{parSigma} to the MLE object.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSSharveyobsFI}} \code{\link{MARSShessian.numerical}} \code{\link{MARSSparamCIs}} \code{\link{marssMLE}} }
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11), ]
MLEobj <- MARSS(dat)
MLEobj.hessian <- MARSShessian(MLEobj)

# show the approx Hessian
MLEobj.hessian$Hessian

# generate a parameter sample using the Hessian
# this uses the rmvnorm function in the mvtnorm package
hess.params <- mvtnorm::rmvnorm(1,
  mean = MLEobj.hessian$parMean,
  sigma = MLEobj.hessian$parSigma
)
}

