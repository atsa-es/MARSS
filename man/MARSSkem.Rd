\name{MARSSkem}
\alias{MARSSkem}
\title{ Maximum Likelihood Estimation for Multivariate Autoregressive State-Space Models }
\description{
  \code{MARSSkem()} performs maximum-likelihood estimation, using an EM algorithm for constrained and unconstrained MARSS models.  Users would not call this function directly normally.  The function \code{\link{MARSS}} calls MARsskem.  However users might want to used MARSSkem directly if they need to avoid some of the error-checking overhead associated with the \code{\link{MARSS}} function.
}
\usage{
MARSSkem(MLEobj)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. }
}
\details{
Objects of class \code{\link{marssMLE}} may be built from scatch but are easier to construct using \code{\link{MARSS}} with \code{MARSS(..., fit=FALSE)}.

  Options for \code{MARSSkem()} may be set using \code{MLEobj$control}.  The commonly used elements of \code{control} are follows (see \code{\link{marssMLE}}:
  \describe{
      \item{\code{minit}}{ Minimum number of EM iterations. You can use this to force the algorithm to do a certain number of iterations.  This is helpful if your soln is not converging.} 
      \item{\code{maxit}}{ Maximum number of EM iterations. } 
      \item{\code{min.iter.conv.test}}{ The minimum number of iterations before the log-log convergence test will be computed.  If \code{maxit} is set less than this, then convergence will not be computed (and the algorithm will just run for maxit iterations). } 
      \item{\code{kf.x0}}{ Whether to set the prior at t=0 ("x00") or at t=1 ("x10"). The default is "x00". } 
      \item{\code{conv.test.deltaT}}{ The number of iterations to use in the log-log convergence test.  This defaults to 9. } 
      \item{\code{abstol}}{ Tolerance for log-likelihood change for the delta logLik convergence test.  If log-likelihood changes less than this amount relative to the previous iteration, the EM algorithm exits. This is normally (default) set to NULL and the log-log convergence test is used instead.} 
      \item{\code{allow.degen}}{ Whether to try setting Q or R elements to zero if they appear to be going to zero. }
      \item{\code{trace}}{ A positive integer. If not 0, a record will be created of each variable over all EM iterations and detailed warning messages (if appropriate) will be printed. }
      \item{\code{safe}}{ If TRUE, \code{MARSSkem} will rerun \code{\link{MARSSkf}} after each individual parameter update rather than only after all parameters are updated.  The latter is slower and unnecessary for many models, but in some cases, the safer and slower algorithm is needed because the ML parameter matrices have high condition numbers. }
      \item{\code{silent}}{ Suppresses printing of progress bars, error messages, warnings and convergence information. }    
    }
}
\value{
  The \code{\link{marssMLE}} object which was passed in, with additional components:
  \item{method}{String "kem".}
  \item{kf}{Kalman filter output. }
  \item{iter.record}{If \code{MLEobj$control$trace = TRUE}, a list with \code{par} = a record of each estimated parameter over all EM iterations and \code{logLik} = a record of the log likelikelihood at each iteration. }
  \item{numIter}{Number of iterations needed for convergence.}
  \item{convergence}{ Did estimation converge successfully? 
    \describe{
      \item{convergence=0}{ Converged in both the abstol test and the log-log plot test. }  
      \item{convergence=1}{Some of the parameter estimates did not converge (based on the log-log plot test AND abstol tests) before \code{MLEobj$control$maxit} was reached. This is not an error per se.}
      \item{convergence=2}{No convergence diagnostics were computed because the MLE object had problems and was not fit. This isn't a convergence error just information.}
      \item{convergence=3}{No convergence diagnostics were computed because the MLE object was not fit. This isn't a convergence error just information.}   
      \item{convergence=10}{Abstol convergence only. Some of the parameter estimates did not converge (based on the log-log plot test) before \code{MLEobj$control$maxit} was reached.  However \code{MLEobj$control$abstol} was reached. }
      \item{convergence=11}{Log-log convergence only. Some of the parameter estimates did not converge (based on the abstol test) before \code{MLEobj$control$maxit} was reached.  However the log-log convergence test was passed. }
      \item{convergence=12}{Abstol convergence only. Log-log convergence test was not computed because \code{MLEobj$control$maxit} was set to less than \code{control$min.iter.conv.test}. }  
      \item{convergence=13}{Lack of convergence info. Parameter estimates did not converge based on the abstol test before \code{MLEobj$control$maxit} was reached.  No log-log information since \code{control$min.iter.conv.test} is less than \code{MLEobj$control$maxit} so no log-log plot test could be done. }
      \item{convergence=42}{\code{MLEobj$control$abstol} was reached but the log-log plot test returned NAs.  This is an odd error and you should set \code{control$trace=TRUE} and look at the outputted \code{$iter.record} to see what is wrong. }
      \item{convergence=52}{The EM algorithm was abandoned due to numerical errors.  Usually this means one of the variances either went to zero or to all elements being equal.  This is not an error per se.  Most likely it means that your model is not very good for your data (too inflexible or too many parameters).  Try setting \code{control$trace=1} to view a detailed error report.}
      \item{convergence=62}{The algorithm was abandoned due to errors in the log-log convergence test.  You should not get this error (it is included for debugging purposes to catch improper arguments passed into the log-log convergence test). }
      \item{convergence=63}{The algorithm was run for \code{control$maxit} iterations, \code{control$abstol} not reached, and the log-log convergence test returned errors.  You should not get this error (it is included for debugging purposes to catch improper arguments passed into the log-log convergence test).}
      \item{convergence=72}{Other convergence errors.  This is included for debugging purposes to catch misc. errors. }        
      }
  }
  \item{logLik}{Log-likelihood.}
  \item{states}{State estimates from the Kalman filter.}
  \item{states.se}{Confidence intervals based on state standard errors, see caption of Fig 6.3 (p. 337) Shumway & Stoffer.}
  \item{errors}{Any error messages.}
}

\section{Discussion}{
 To ensure that the global maximum-likelihood values are found, it is recommended that you test the fit under different initial parameter values, particularly if the model is not a good fit to the data.  This requires more compuation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.   For many models and for draft analyses, this is unnecessary, but answers should be checked using an initial conditions search before reporting final values.  See the chapter on initial conditions in the User Guide for a discussion on how to do this.

 \code{MARSSkem()} calls a Kalman filter/smoother (\code{\link{MARSSkf}}) for hidden state estimation.   The algorithm allows two options for the initial state conditions: fixed but unknown or a prior.  In the first case, x0 (whether at t=0 or t=1) is treated as fixed but unknown (estimated); in this case, \code{fixed$V0=0} and x0 is estimated.  This is the default behavior.  In the second case, the initial conditions are specified with a prior and V0!=0.   In the later case, x0 or V0 may be estimated.  MARSS will allow you to try to estimate both, but many researchers have noted that this is not robust so you should fix one or the other.
 
 If you get errors, it generally means that the solution involves an ill-conditioned matrix.    For example, your Q or R matrix is going to a value in which all elements have the same value, for example zero.  If for example, you tried to fit a model with fixed and high R matrix and the variance in that R matrix was much higher than what is actually in the data, then you might drive Q to zero.   Also if you try to fit a structurally inadequate model, then it is not unusual that Q will be driven to zero.  For example, if you fit a model with 1 hidden state trajectory to data that clearly have 2 quite different hidden state trajectories, you might have this problem.  Comparing the likelihood of this model to a model with more structural flexibility should reveal that the structually inflexible model is inadequate (much lower likelihood).
 
 Convergence testing is done via a combination of two tests.  The first test (abstol test) is the test that the change in the absolute value of the log-likelihood from one iteration to another is less than some tolerance value (abstol).  The second test (log-log test) is that the slope of a plot of the log of the parameter value or log-likelihood versus the log of the iteration number is less than some tolerance.  Both of these must be met to generate the Success! parameters converged output.  If you want to circumvent one of these tests, then set the tolerance for the unwanted test to be high.  That will guarantee that that test is met before the convergence test you want to use is met.  The tolerance for the abstol test is set by \code{control$abstol} and the tolerance for the log-log test is set by \code{control$conv.test.slope.tol}.  Anything over 1 is huge for both of these. 
     }

\references{ 
  R. H. Shumway and D. S. Stoffer (2006).  Chapter 6 in Time series analysis and its applications.  Springer-Verlag, New York.
  
  Ghahramani, Z. and Hinton, G. E. (1996) Parameter estimation for linear dynamical systems. Technical Report CRG-TR-96-2, University of
Totronto, Dept. of Computer Science. 

Harvey, A. C. (1989) Chapter 5 in Forecasting, structural time series models and the
Kalman filter. Cambridge University Press, Cambridge, UK.
  
The user guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.      

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive
state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]\code{RShowDoc("EMDerivation",package="MARSS")} to open a copy.
}
\author{ 
  Eli Holmes and Eric Ward, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov
}
\seealso{ 
  \code{\link{MARSSkf}}, \code{\link{marssMLE}}, \code{\link{MARSSoptim}}
  }  
\examples{
  dat = t(harborSeal)
  dat = dat[2:4,] 
  #you can use MARSS to construct a proper marssMLE object.
  MLEobj = MARSS(dat, model=list(Q="diagonal and equal", U="equal"), fit=FALSE) 
  #Pass this MLEobj to MARSSkem to do the fit.
  kemfit = MARSSkem(MLEobj)
}
