\name{MARSS-package}
\alias{MARSS-package}
\docType{package}
\title{
Multivariate Autoregressive State-Space Model Estimation
}
\description{
The MARSS package fits constrained and unconstrained multivariate autoregressive time-series models to multivariate time series data.  To open the user guide from the command line, type \code{RShowDoc("UserGuide",package="MARSS")}.  To open an overview page with package information and links to the R scripts in the User Guide, type \code{RShowDoc("index",package="MARSS")}.    

The main function is \code{\link{MARSS}} which is used to fit a specified model to data and estimate the model parameters.  MARSS model specification is based on "form" (an argument to a MARSS() call).  The form tells MARSS() what to expect in the model list (model is a MARSS argument) and how to translate that into the base model form used in the fitting algorithms.

The default MARSS model form is "marxss", which is a model of the following form:
\describe{
\item{x(t+1) = B x(t) + U + C c(t) + G w(t), where w(t) ~ MVN(0,Q)}{}
\item{y(t) = Z x(t) + A + D d(t) + H v(t), where v(t) ~ MVN(0,R)}{}
\item{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0) }{} 
}
The parameters, hidden state processes (x), and observations (y) are matrices:
  \itemize{
    \item x(t) is m x 1 
    \item y(t) is n x 1 (m<=n)
    \item Z is n x m 
    \item B is m x m 
    \item U is m x 1 
    \item Q is g x g (typically m x m) 
    \item G is m x g
    \item A is n x 1 
    \item R is h x h (typically n x n)
    \item H is n x h
    \item C is m x q
    \item D is n x p
    \item c(t) is q x 1
    \item d(t) is q x 1
    \item x0 is m x 1 
    \item V0 is m x m 
  }
  All parameters can be time-varying.
  
The package functions estimate the model parameters using an EM algorithm (primarily but see \code{\link{MARSSoptim}}).  Parameters may be constrained to have shared elements (elements which are constrained to have the same value) or fixed elements (with the other elements estimated).  The states and smoothed state estimates are provided via a Kalman filter and smoother.  Bootstrapping, confidence interval estimation, bias estimation, model selection and simulation functions are provided.   The main user interface to the package is the top-level function \code{\link{MARSS}}.  
}
\details{
Important MARSS functions:
\describe{
      \item{\code{\link{MARSS}}}{ Top-level function for specifying and fitting MARSS models. } 
      \item{\code{\link{MARSSsimulate}}}{ Produces simulated data from a MARSS model. } 
      \item{\code{\link{MARSSkem}}}{ Estimates MARSS parameters using an EM algorithm. } 
      \item{\code{\link{MARSSkf}} }{ Kalman filter and smoother. } 
      \item{\code{\link{MARSSoptim}}}{ Estimates MARSS parameters using a quasi-Newton algorithm via \code{\link{optim}}. } 
      \item{\code{\link{MARSSaic}}}{ Calculates AICc, AICc, and various bootstrap AICs. } 
      \item{\code{\link{MARSSboot}}}{ Creates bootstrap MARSS parameter estimates. }
      \item{\code{\link{MARSSparamCIs}}}{ Computes confidence intervals for maximum-likelihood estimates of MARSS parameters. }    
    }
}
\author{
  Eli Holmes, Eric Ward and Kellie Wills, NOAA, Seattle, USA.

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov, 
  
  kellie(dot)wills(at)noaa(dot)gov
}
\references{
The MARSS user guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.

Type \code{RShowDoc("index",package="MARSS")} to see all the package documentation, tutorials, and R scripts from the User Guide.
}

\keyword{ package }
