\name{MARSSoptim}
\alias{MARSSoptim}
\title{Parameter estimation for MARSS models using optim}
\description{
  Parameter estimation for MARSS models using R's \code{\link{optim}} function.   This allows access to R's quasi-Newton algorithms available via the \code{\link{optim}} function.  The \code{MARSSoptim} function is called when \code{\link{MARSS}} is called with \code{method="BFGS"}.  This is a base function in the \code{\link{MARSS-package}}. 
}
\usage{
MARSSoptim(MLEobj)
}
\arguments{
  \item{MLEobj}{An object of class \code{\link{marssMLE}}.}
}
\details{
Objects of class \code{\link{marssMLE}} may be built from scratch but are easier to construct using \code{\link{MARSS}} with \code{MARSS(..., fit=FALSE, method="BFGS")}.

  Options for \code{\link{optim}} are passed in using \code{MLEobj$control}.  See \code{\link{optim}} for a list of that function's control options.  If \code{lower} and \code{upper} for \code{\link{optim}} need to be passed in, they should be passed in as part of \code{control} as \code{control$lower} and \code{control$upper}. Additional \code{control} arguments affect printing and initial conditions.
  \describe{
      \item{\code{MLEobj$control$kf.x0}}{ The initial condition is at $t=0$ if kf.x0="x00".  The initial condition is at $t=1$ if kf.x0="x10". }    
      \item{\code{MLEobj$marss$diffuse}}{ If diffuse=TRUE, a diffuse initial condition is used. MLEobj$par$V0 is then the scaling function for the diffuse part of the prior.  Thus the prior is V0*kappa where kappa-->Inf.  Note that setting a diffuse prior does not change the correlation structure within the prior. If diffuse=FALSE, a non-diffuse prior is used and MLEobj$par$V0 is the non-diffuse prior variance on the initial states. The the prior is V0.}      \item{\code{MLEobj$control$silent}}{ Suppresses printing of progress bars, error messages, warnings and convergence information. }    
    }
}
\value{
  The \code{\link{marssMLE}} object which was passed in, with additional components:
  \item{method}{String "BFGS".}
  \item{kf}{Kalman filter output. }
  \item{iter.record}{If \code{MLEobj$control$trace = TRUE}, then this is the \code{$message} value from \code{\link{optim}}.}
  \item{numIter}{Number of iterations needed for convergence.}
  \item{convergence}{ Did estimation converge successfully? 
    \describe{
      \item{convergence=0}{ Converged in less than \code{MLEobj$control$maxit} iterations and no evidence of degenerate solution. } 
      \item{convergence=1}{ Maximum number of iterations \code{MLEobj$control$maxit} was reached before \code{MLEobj$control$abstol} condition was satisfied. } 
      \item{convergence=10}{ Some of the variance elements appear to be degenerate.  T} 
      \item{convergence=52}{ The algorithm was abandoned due to errors from the "L-BFGS-B" method.  } 
      \item{convergence=53}{ The algorithm was abandoned due to numerical errors in the likelihood calculation from \code{\link{MARSSkf}}.  If this happens with "BFGS", it can sometimes be helped with a better initial condition.  Try using the EM algorithm first (\code{method="kem"}), and then using the parameter estimates from that to as initial conditions for \code{method="BFGS"}.  } 
      }
  }
  \item{logLik}{Log-likelihood.}
  \item{states}{State estimates from the Kalman filter.}
  \item{states.se}{Confidence intervals based on state standard errors, see caption of Fig 6.3 (p. 337) Shumway & Stoffer.}
  \item{errors}{Any error messages.}
}
\section{Discussion}{
 The function only returns parameter estimates.  To compute CIs, use \code{\link{MARSSparamCIs}} but if you use parametric or non-parametric bootstrapping with this function, it will use the EM algorithm to compute the bootstrap parameter estimates!  The quasi-Newton estimates are too fragile for the bootstrap routine since one often needs to search to find a set of initial conditions that work (i.e. don't lead to numerical errors).
 
 Estimates from \code{MARSSoptim} (which come from \code{\link{optim}}) should be checked against estimates from the EM algorithm.  If the quasi-Newton algorithm works, it will tend to find parameters with higher likelihood faster than the EM algorithm.  However, the MARSS likelihood surface can be multimodal with sharp peaks at degenerate solutions where a Q or R diagonal element equals 0.  The quasi-Newton algorithm sometimes gets stuck on these peaks even when they are not the maximum.  Neither an initial conditions search nor starting near the known maximum (or from the parameters estimates after the EM algorithm) will necessarily solve this problem.   Thus it is wise to check against EM estimates to ensure that the BFGS estimates are close to the MLE estimates (and vis-a-versa, it's wise to rerun with method="BFGS" after using method="kem").  Conversely, there is a strong flat ridge in your likelihood, the EM algorithm can report early convergence while the BFGS may continue much further along the ridge and find very different parameter values.  Of course a likelihood surface with strong flat ridges makes the MLEs less informative...
 
 Note this is mainly a problem if the time series are short or very gappy.  If the time series are long, then the likelihood surface should be nice with a single interior peak.  In this case, the quasi-Newton algorithm works well but it can still be sensitive (and slow) if not started with a good initial condition.  Thus starting it with the estimates from the EM algorithm is often desirable.
 
 One should be aware that the prior set on the variance of the initial states at t=0  or t=1 can have catastrophic effects on one's estimates if the presumed prior covariance structure conflicts with the structure implied by the MARSS model.  For example, if you use a diagonal variance-covariance matrix for the prior but the model implies a matrix with non-zero covariances, your MLE estimates can be strongly influenced by the prior variance-covariance matrix.  Setting a diffuse prior does not help because the diffuse prior still has the correlation structure specified by V0.  One way to detect priors effects is to compare the BFGS estimates to the EM estimates.  Persistent differences typically signify a problem with the correlation structure in the prior conflicting with the implied correlation structure in the MARSS model.  If this is the case, using V0=0 and estimating the x0 elements (with control$kf.x0="x10") can often help.
}
\author{
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{  \code{\link{MARSS}} \code{\link{MARSSkem}}  \code{\link{marssMLE}}  \code{\link{optim}}
}
\examples{
dat <- t(harborSealWA)
dat <- dat[2:4, ] # remove the year row

# fit a model with EM and then use that fit as the start for BFGS
# fit a model with 1 hidden state where obs errors are iid
# R="diagonal and equal" is the default so not specified
# Q is fixed
kemfit <- MARSS(dat, model = list(Z = matrix(1, 3, 1), Q = matrix(.01)))
bfgsfit <- MARSS(dat,
  model = list(Z = matrix(1, 3, 1), Q = matrix(.01)),
  inits = coef(kemfit, form = "marss"), method = "BFGS"
)
}
