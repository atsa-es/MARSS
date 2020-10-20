\name{MARSShessian}
\alias{MARSShessian}

\title{ Parameter Variance-Covariance Matrix from the Hessian Matrix}
\description{
Calculates an approximate parameter variance-covariance matrix for the parameters using an inverse of the Hessian of the negative log-likelihood function at the MLEs (the observed Fisher Information matrix).   It appends \code{$Hessian}, \code{$parMean}, \code{$parSigma} to the \code{\link{marssMLE}} object. 
}
\usage{
MARSShessian(MLEobj, method=c("Harvey1989", "fdHess", "optim"))
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. 
  This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{\link{MARSSkem}}. }
  \item{method}{ The method to use for computing the Hessian. Options are \code{Harvey1989} to use the Harvey (1989) recursion, which is an analytical solution, \code{fdHess} or \code{optim} which are two numerical methods.  Although \code{optim} can be passed to this function, in the internal functions which call this function, \code{fdHess} will be used if a numerical estimate is requested.}
 }
\details{
See \code{\link{MARSSFisherI}} for a discussion of the observed Fisher Information matrix and references.

   Method \code{fdHess} uses \code{\link[nlme]{fdHess}} from package nlme to numerically estimate the Hessian matrix (the matrix of partial 2nd derivatives of the negative log-likelihood function at the MLE).     Method \code{optim} uses \code{\link{optim}} with \code{hessian=TRUE} and \code{list(maxit=0)} to ensure that the Hessian is computed at the values in the \code{par} element of the MLE object. Method \code{Harvey1989} (the default) uses the recursion in Harvey (1989) to compute the observed Fisher Information of a MARSS model analytically.
   
   Note that the parameter confidence intervals computed with the observed Fisher Information matrix are based on the asymptotic normality of maximum-likelihood estimates under a large-sample approximation. 
   
}
\value{
  \code{MARSShessian()} attaches 
   \code{Hessian}, \code{parMean} and \code{parSigma} to the \code{\link{marssMLE}} object that is passed into the function.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
} 
\seealso{ \code{\link{MARSSFisherI}}, \code{\link{MARSSharveyobsFI}}, \code{\link{MARSShessian.numerical}}, \code{\link{MARSSparamCIs}}, \code{\link{marssMLE}} }
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
\references{ 
Harvey, A. C. (1989) Section 3.4.5 (Information matrix) in Forecasting, structural time series models and the Kalman filter. Cambridge University Press, Cambridge, UK.

See also J. E. Cavanaugh and R. H. Shumway (1996) On computing the expected Fisher information matrix for state-space model parameters. Statistics & Probability Letters 26: 347-355.  This paper discusses the Harvey (1989) recursion (and proposes an alternative).
}

