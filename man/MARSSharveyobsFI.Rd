\name{MARSSharveyobsFI}
\alias{MARSSharveyobsFI}
\title{ MARSS Hessian Matrix via the Harvey (1989) Recursion}
\description{
Calculates the observed Fisher Information analytically via the recursion by Harvey (1989).  This is the same as the Hessian of the log-likelihood function at the MLEs.  This is a utility function in the \code{\link{MARSS-package}} and is not exported.  Use \code{\link{MARSShessian}} to access.
}
\usage{
MARSSharveyobsFI( MLEobj )
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. 
  This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{\link{MARSSkem}}. }
 }
\value{
  The observed Fisher Information matrix computed via equation 3.4.69 in Harvey (1989).  The differentials in the equation are computed in the recursion in equations 3.4.73a to 3.4.74b.  Harvey (1989) discusses missing observations in section 3.4.7.  However, the MARSSharveyobsFI() function implements the approach of Shumway and Stoffer (2006) in section 6.4 for the missing values.  See Holmes (2012) for a full discussion of the missing values modifications.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSShessian}} \code{\link{MARSSparamCIs}} }
\examples{
  dat = t(harborSeal)
  dat = dat[c(2,11),]
  MLEobj = MARSS(dat)
  MARSS:::MARSSharveyobsFI(MLEobj)
}
\references{ 
R. H. Shumway and D. S. Stoffer (2006).  Section 6.4 (Missing Data Modifications) in Time series analysis and its applications.  Springer-Verlag, New York.
  
Harvey, A. C. (1989) Section 3.4.5 (Information matrix) in Forecasting, structural time series models and the Kalman filter. Cambridge University Press, Cambridge, UK.
  
Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive
state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]\code{RShowDoc("EMDerivation",package="MARSS")} to open a copy.

See also J. E. Cavanaugh and R. H. Shumway (1996) On computing the expected Fisher information matrix for state-space model parameters. Statistics & Probability Letters 26: 347-355.  This paper discusses the Harvey (1989) recursion (and proposes an alternative).
}