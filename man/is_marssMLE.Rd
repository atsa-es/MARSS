\name{is.marssMLE}
\alias{is.marssMLE}

\title{ Tests marssMLE object for completeness }
\description{
  Tests an marssMLE for completeness to determine if it has all the pieces and attributes necessary to be passed to MARSS functions for fitting, filtering, smoothing, or displaying. Internal function, use \code{MARSS:::} to access. This is a very slow function which should not be called repeatedly in a \code{for} loop for example.
}
\usage{
is.marssMLE(MLEobj)
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. See Details.}
}
\details{
  The \code{is.marssMLE()} function checks components \code{marss}, \code{start} and \code{control}, which must be present for estimation by functions e.g. \code{\link{MARSSkem}()}. Components returned from estimation must include at least \code{method}, \code{par}, \code{kf}, \code{numIter}, \code{convergence} and \code{logLik}. Additional components (e.g. AIC) may be returned, as described in function help files.
\describe{
  \item{\code{model}}{ An object of class \code{\link{marssMODEL}} in whatever form the user specified in the call to \code{\link{MARSS}()}.  Default is form "marxss". }
  \item{\code{marss}}{ An object of class \code{\link{marssMODEL}} in "marss" forms, needed for all the base MARSS functions. }
  \item{\code{start}}{ List with matrices specifying initial values for parameters to be used (if needed) by the maximization algorithm. 
    \describe{
    \item{\code{B}}{ Initial value(s) for B matrix (m x m). }
    \item{\code{U}}{ Initial value(s) for U matrix (m x 1). }
    \item{\code{Q}}{ Initial value(s) for Q variance-covariance matrix (m x m). }
    \item{\code{Z}}{ Initial value(s) for Z matrix (n x m). }
    \item{\code{A}}{ Initial value(s) for A matrix (n x 1). }
    \item{\code{R}}{ Initial value(s) for R variance-covariance matrix (n x n). }
    \item{\code{x0}}{ Initial value(s) for initial state vector (m x 1). }
    \item{\code{V0}}{ Initial variance(s) for initial state variance (m x m). }
    }
  }
  \item{\code{control}}{ A list specifying estimation options.  The following options are needed by \code{\link{MARSSkem}()}.  Other control options
  can be set if needed for other estimation methods, e.g. the control options listed for \code{\link{optim}} for use with \code{\link{MARSSoptim}()}.  The default values for control options are set in \code{alldefaults[[method]]} which is specified in \code{MARSSsettings.R}.
    \describe{
      \item{\code{minit}}{ The minimum number of iterations to do in the maximization routine (if needed by method). } 
      \item{\code{maxit}}{ Maximum number of iterations to be used in the maximization routine (if needed by method). }
      \item{\code{min.iter.conv.test}}{ Minimum iterations to run before testing convergence via the slope of the log parameter versus log iterations.} 
      \item{\code{conv.test.deltaT=9}}{ Number of iterations to use for the testing convergence via the slope of the log parameter versus log iterations. }
      \item{\code{conv.test.slope.tol}}{ The slope of the log parameter versus log iteration to use as the cut-off for convergence.  The default is 0.5 which is a bit high.  For final analyses, this should be set lower. }
      \item{\code{abstol}}{ The logLik.(iter-1)-logLik.(iter) convergence tolerance for the maximization routine. Both the abstol and the slope of the log of the parameters versus the log iteration tests must be met for convergence.} 
      \item{\code{trace}}{ A positive integer.  If not 0, a record will be created during maximization iterations (what's recorded depends on method of maximization). -1 turns off most internal error checking.}
      \item{\code{safe}}{ Logical.  If TRUE, then the Kalman filter is run after each update equation in the EM algorithm.  This slows down the algorithm.  The default is FALSE. }
      \item{\code{allow.degen}}{ If TRUE, replace Q or R diagonal elements by 0 when they become very small.} 
      \item{\code{min.degen.iter}}{ Number of iterations before trying to set a diagonal element of Q or R to zero). } 
      \item{\code{degen.lim}}{ How small the Q or R diagonal element should be before attempting to replace it with zero. }
  \item{\code{silent}}{ Suppresses printing of progress bar, error messages and convergence information. } 
    }
  }
  \item{\code{method}}{ A string specifying the estimation method. MARSS allows "kem" for EM and "BFGS" for quasi-Newton. }
Once the model has been fitted, additional elements are added.
  \item{\code{par}}{ A list with 8 matrices of estimated parameter values Z, A, R, B, U, Q, x0, V0. }
  \item{\code{states}}{ Expected values of the x (hidden states). }
  \item{\code{states.se}}{ Standard errors on the estimates states. }
  \item{\code{ytT}}{ Expected values of the y.  This is just y for non-missing y. }
  \item{\code{ytT.se}}{ Standard errors on the ytT. This will be 0 for non-missing y. }
  \item{\code{kf}}{ A list containing Kalman filter/smoother output if \code{control$trace} is > 0.}
  \item{\code{Ey}}{ A list containing expectations involving y. Output if \code{control$trace} is > 0.}
  \item{\code{numIter}}{ Number of iterations which were required for convergence. }
  \item{\code{convergence}}{ Convergence status and errors. 0 means converged successfully.  Anything else means an error or warning.}
  \item{\code{logLik}}{ Log-likelihood. }
  \item{\code{AIC}}{ AIC }
  \item{\code{AICc}}{ Corrected AIC. }
  \item{\code{call}}{ A list of all the arguments passed into the MARSS call.  Not required for most functions, but is a record of what was used to call MARSS for checking and can be used to customize the printing of MARSS output. }
}
}

\value{
  TRUE if no problems; otherwise a message describing the problems.
}

\author{ 
  Kellie Wills, NOAA, Seattle, USA.  

  kellie(dot)wills(at)noaa(dot)gov
}
\seealso{ 
  \code{\link{marssMODEL}}, \code{\link{MARSSkem}} 
}


