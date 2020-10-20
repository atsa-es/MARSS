\name{MARSSFisherI}
\alias{MARSSFisherI}

\title{ Observed Fisher Information Matrix at the MLE }
\description{
Returns the observed Fisher Information matrix for a \code{\link{marssMLE}} object (a fitted MARSS model) via either the analytical algorithm of Harvey (1989) or a numerical estimate. 

The observed Fisher Information is the negative of the second-order partial derivatives of the log-likelihood function evaluated at the MLE. The derivatives being with respect to the parameters. The Hessian matrix is the second-order partial derivatives of a scalar-valued function. Thus the observed Fisher Information matrix is the Hessian of the negative log-likelihood function evaluated at the MLE (or equivalently the negative of the Hessian of the log-likelihood function). The inverse of the observed Fisher Information matrix is an estimate of the asymptotic variance-covariance matrix for the estimated parameters. Use \code{\link{MARSShessian}()} (which calls \code{MARSSFisherI()}) to return the parameter variance-covariance matrix computed from the observed Fisher Information matrix.

Note for the numerically estimated Hessian, we pass in the negative log-likelihood function to a minimization function. As a result, the numerical functions return the Hessian of the negative log-likelihood function (which is the observed Fisher Information matrix).
}
\usage{
MARSSFisherI(MLEobj, method=c("Harvey1989", "fdHess", "optim"))
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. 
  This object must have a \code{$par} element containing MLE parameter estimates from e.g. \code{\link{MARSSkem}()}. }
  \item{method}{ The method to use for computing the observed Fisher Information matrix. Options are \code{"Harvey1989"} to use the Harvey (1989) recursion, which is an analytical solution, \code{"fdHess"} or \code{"optim"} which are two numerical methods.  Although 'optim' can be passed to the function, 'fdHess' is used for all numerical estimates used in the MARSS package.}
 }
\details{
   Method 'fdHess' uses \code{\link{fdHess}()} from package \code{\link{nlme}} to numerically estimate the Hessian of the negative log-likelihood function at the MLEs. Method 'optim' uses \code{\link{optim}()} with \code{hessian=TRUE} and \code{list(maxit=0)} to ensure that the Hessian is computed at the values in the \code{par} element of the MLE object. The \code{par} element of the \code{\link{marssMLE}} object is the MLE.  
   
Method 'Harvey1989' (the default) uses the recursion in Harvey (1989) to compute the observed Fisher Information of a MARSS model analytically. See Holmes (2016c) for a discussion of the Harvey (1989) algorithm and see Holmes (2017) on how to implement the algorithm for MARSS models with linear constraints (the type of MARSS models that the MARSS R package addresses).
   
There has been research on computing the observed Fisher Information matrix from the derivatives used by EM algorithms (discussed in Holmes (2016a, 2016b)), for example Louis (1982). Unfortunately, the EM algorithm used in the MARSS package is for time series data and the temporal correlation must be dealt with, e.g. Duan & Fulop (2011). Oakes (1999) has an approach that only involves derivatives of \eqn{E[LL(\Theta)|\mathbf{y},\Theta']}{E(LL(Theta)|data, Theta')} but one of the derivatives will be the derivative of the \eqn{E[\mathbf{X}|\mathbf{y},\Theta']}{E(X|data, Theta')} with respect to \eqn{\Theta'}{Theta'}.  It is not clear how to do that derivative.  Moon-Ho, Shumway and Ombao (2006) suggest (page 157) that this derivative is hard to compute.

}
\value{
  Returns the observed Fisher Information matrix.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
} 
\seealso{ \code{\link{MARSSharveyobsFI}}, \code{\link{MARSShessian.numerical}}, \code{\link{MARSSparamCIs}}, \code{\link{marssMLE}} }
\examples{
dat <- t(harborSeal)
dat <- dat[2:4, ]
MLEobj <- MARSS(dat, model=list(Z=matrix(1,3,1), R="diagonal and equal"))
MARSSFisherI(MLEobj)
MARSSFisherI(MLEobj, method="fdHess")

}
\references{ 
Harvey, A. C. (1989) Section 3.4.5 (Information matrix) in Forecasting, structural time series models and the Kalman filter. Cambridge University Press, Cambridge, UK.
  
See also J. E. Cavanaugh and R. H. Shumway (1996) On computing the expected Fisher information matrix for state-space model parameters. Statistics & Probability Letters 26: 347-355.  This paper discusses the Harvey (1989) recursion (and proposes an alternative).

Holmes, E. E. 2016a. Notes on computing the Fisher Information matrix for MARSS models. Part I Background. Technical Report. https://doi.org/10.13140/RG.2.2.27306.11204/1  \href{https://eeholmes.github.io/posts/2016-5-18-FI-recursion-1/}{Notes}

Holmes, E. E. 2016b. Notes on computing the Fisher Information matrix for MARSS models. Part II Louis 1982. Technical Report. https://doi.org/10.13140/RG.2.2.35694.72000  \href{https://eeholmes.github.io/posts/2016-5-19-FI-recursion-2/}{Notes}

Holmes, E. E. 2016c. Notes on computing the Fisher Information matrix for MARSS models. Part III Overview of Harvey 1989. https://eeholmes.github.io/posts/2016-6-16-FI-recursion-3/

Holmes, E. E. 2017. Notes on computing the Fisher Information matrix for MARSS models. Part IV Implementing the Recursion in Harvey 1989. https://eeholmes.github.io/posts/2017-5-31-FI-recursion-4/

Duan, J. C. and A. Fulop. (2011) A stable estimator of the information matrix under EM for dependent data. Statistics and Computing 21: 83-91

Louis, T. A. 1982. Finding the observed information matrix when using the EM algorithm. Journal of the Royal Statistical Society. Series B (Methodological). 44: 226-233.

Oakes, D. 1999. Direct calculation of the information matrix via the EM algorithm. Journal of the Royal Statistical Society. Series B (Methodological). 61: 479-482.

Moon-Ho, R. H., R. H. Shumway, and Ombao 2006. The state-space approach to modeling dynamic processes. Chapter 7 in Models for Intensive Longitudinal Data. Oxford University Press.
}
