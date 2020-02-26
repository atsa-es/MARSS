\name{MARSSparamCIs}
\alias{MARSSparamCIs}
\title{ Standard Errors, Confidence Intervals and Bias for MARSS Parameters }
\description{
Computes standard errors, confidence intervals and bias for the maximum-likelihood estimates of MARSS model parameters.  If you want confidence intervals on the estimated hidden states, see \code{\link{print.marssMLE}} and look for "states.cis".
}
\usage{
MARSSparamCIs(MLEobj, method = "hessian", alpha = 0.05, nboot =
    1000, silent = TRUE, hessian.fun = "Harvey1989")
}

\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. Must have a \code{$par} element containing the MLE parameter estimates. }
  \item{method}{ Method for calculating the standard errors: "hessian", "parametric", and "innovations" implemented currently. }
  \item{alpha}{ alpha level for the 1-alpha confidence intervals. }
  \item{nboot}{ Number of bootstraps to use for "parametric" and "innovations" methods. }
  \item{hessian.fun}{ The function to use for computing the Hessian.  See \code{\link{MARSShessian}}.}
  \item{silent}{ If false, a progress bar is shown for "parametric" and "innovations" methods. }
}
\details{
  Approximate confidence intervals (CIs) on the model parameters may be calculated from the observed Fisher Information matrix ("Hessian CIs", see \code{\link{MARSSFisherI}}) or parametric or non-parametric (innovations) bootstrapping using \code{nboot} bootstraps.  The Hessian CIs are based on the asymptotic normality of MLE parameters under a large-sample approximation.  The Hessian computation for variance-covariance matrices is a symmetric approximation and the lower CIs for variances might be negative.  Bootstrap estimates of parameter bias are reported if method "parametric" or "innovations" is specified.
  
Note, these are added to the \code{par}  elements of a \code{\link{marssMLE}} object but are in "marss" form not "marxss" form. Thus the \code{MLEobj$par.upCI} and related elements that are added to the \code{\link{marssMLE}} object may not look familiar to the user.  Instead the user should extract these elements using \code{print(MLEobj)} and passing in the argument \code{what} set to "par.se","par.bias","par.lowCIs", or "par.upCIs".  See \code{\link{print.marssMLE}}. Or use \code{\link{tidy.marssMLE}}.
}
\value{
  \code{MARSSparamCIs} returns the \code{\link{marssMLE}} object passed in, with additional components \code{par.se}, \code{par.upCI}, \code{par.lowCI}, \code{par.CI.alpha}, \code{par.CI.method}, \code{par.CI.nboot} and \code{par.bias} (if method is "parametric" or "innovations").
}
\references{ 
Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.
 }
\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ \code{\link{MARSSboot}}, \code{\link{MARSSinnovationsboot}}, \code{\link{MARSShessian}}
 }
\examples{
dat <- t(harborSealWA)
dat <- dat[2:4, ]
kem <- MARSS(dat, model = list(
  Z = matrix(1, 3, 1),
  R = "diagonal and unequal"
))
kem.with.CIs.from.hessian <- MARSSparamCIs(kem)
kem.with.CIs.from.hessian
}

