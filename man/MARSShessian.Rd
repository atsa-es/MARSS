\name{MARSShessian}
\alias{MARSShessian}
\alias{MARSShessian.backtrans}
\title{ MARSS Parameter Variance-Covariance Matrix from the Hessian Matrix}
\description{
Calculates an approximate parameter variance-covariance matrix for the parameters.  The variance-covariance parameters are TRANSFORMED via a Cholesky decomposition so that the variance-covariance matrices remain positive definite.  The function returns a 
  \code{\link{marssMLE}} object in this TRANSFORMED form.  It appends \code{$Hessian}, \code{$parMean}, \code{$parSigma} in this transformed state. This is a utility function in the \code{\link{MARSS-package}}.
}
\usage{
MARSShessian(MLEobj)
MARSShessian.backtrans(MLEobj.hessian, par.hessian)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. 
  This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{\link{MARSSkem}}. }
  \item{MLEobj.hessian}{ An object of class \code{\link{marssMLE}} but with any variance-covariance matrices in their Cholesky decomposition form.  This requires a \code{marssMODEL} with transformed \code{fixed} elements and transformed \code{par} element in the \code{marssMLE} object.}
   \item{par.hessian}{ A vector of parameter estimates where any variance matrix elements are in the Cholesky decomposed form.}
 }
\details{
  Uses \code{\link{fdHess}} from package \code{\link{nlme}} to numerically estimate the Hessian matrix
   (the matrix of partial 2nd derivatives of the parameter estimates).  Hessian CIs are based on 
   the asymptotic normality of ML estimates under a large-sample approximation. 
}
\value{
  \code{MARSShessian()} returns a transformed \code{\link{marssMLE}} object passed in along with additional components
   \code{Hessian}, \code{gradient}, \code{parMean} and \code{parSigma} computed by the \code{MARSShessian} function.
   
  \code{MARSShessian.backtrans(MLEobj.hessian, par.hessian)} takes a transformed \code{\link{marssMLE}} object (as output by \code{MARSShessian} along with a parameter vector in transformed form and returns a parameter vector in regular (non-transformed form).
 
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSSparamCIs}} \code{\link{marssMLE}} }
\examples{
  dat = t(harborSeal)
  dat = dat[c(2,11),]
  MLEobj = MARSS(dat)
  #variance terms are TRANSFORMED by the Cholesky decomposition
  MLEobj.hessian = MARSShessian(MLEobj)

  #show the approx Hessian, again variance terms are transformed
  MLEobj.hessian$Hessian
  
  #generate a (Cholesky transformed) parameter sample using the Hessian
  #this uses the rmvnorm function in the mvtnorm package
  hess.params = mvtnorm::rmvnorm(1, mean=MLEobj.hessian$parMean, 
                           sigma=MLEobj.hessian$parSigma)
        
  #back transform so variance terms are non-transformed
  MARSShessian.backtrans(MLEobj.hessian, hess.params)
}

