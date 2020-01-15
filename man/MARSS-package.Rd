\name{MARSS-package}
\alias{MARSS-package}
\docType{package}
\title{
Multivariate Autoregressive State-Space Model Estimation
}
\description{
The MARSS package fits time-varying constrained and unconstrained multivariate autoregressive time-series models to multivariate time series data.  Get started quickly, go to the \href{https://cran.r-project.org/package=MARSS/vignettes/Quick_Start.pdf}{Quick Start Guide} (or off-line type \code{RShowDoc("Quick_Start",package="MARSS")}). To open the MARSS User Guide with many vignettes and examples, go to \href{https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf}{User Guide} (or off-line, type \code{RShowDoc("UserGuide",package="MARSS")}).  

The default MARSS model form is a MARXSS model: Multivariate Auto-Regressive(1) eXogenous inputs State-Space model. This model has the following form:
\deqn{\mathbf{x}(t+1) = \mathbf{B} x(t) + \mathbf{u} + \mathbf{C} \mathbf{c}(t) + \mathbf{G} \mathbf{w}(t), \textrm{ where } \mathbf{w}(t) \sim \textrm{MVN}(0,\mathbf{Q})}{x(t+1) = B x(t) + u + C c(t) + G w(t), where w(t) ~ MVN(0,Q)}
\deqn{\mathbf{y}(t) = \mathbf{Z} \mathbf{x}(t) + \mathbf{a} + \mathbf{D} \mathbf{d}(t) + \mathbf{H} \mathbf{v}(t), \textrm{ where } \mathbf{v}(t) \sim \textrm{MVN}(0,\mathbf{R})}{y(t) = Z x(t) + a + D d(t) + H v(t), where v(t) ~ MVN(0,R)}
\deqn{\mathbf{X}(1) \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) \textrm{ or } \mathbf{X}(0) \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) }{x(1) ~ \textrm{MVN}(x0, V0) or x(0) ~ \textrm{MVN}(x0, V0) }
All parameters can be time-varying; (t) is left off the parameters to remove clutter. Note, by default \eqn{V0} is a matrix of all zeros and thus \eqn{x(1)} or \eqn{x(0)} is treated as an estimated parameter not a diffuse prior.

The parameter matrices can have fixed values and linear constraints. This is an example of a 3x3 matrix with fixed values and linear constraints. In this example all the matrix elements can be written as a linear function of \eqn{a}, \eqn{b}, and \eqn{c}:
\deqn{\left[\begin{array}{c c c} a+2b & 1 & a\\ 1+3a+b & 0 & b \\ 0 & -2 & c\end{array}\right]}{ [a+2b  1   a \\ 1+3a+b   0   b \\ 0   -2    c ]}
Values such as \eqn{a b} or \eqn{a^2} or \eqn{log(a)} are not allowed as those would not be linear.

The MARSS model parameters, hidden state processes (x), and observations (y) are matrices:
  \itemize{
    \item x(t), x0, and U are m x 1 
    \item y(t) and A are n x 1 (m<=n)
    \item B and V0 are m x m 
    \item Z is n x m 
    \item Q is g x g (default m x m) 
    \item G is m x g (default m x m identity matrix)
    \item R is h x h (default n x n)
    \item H is n x h (default n x n identity matrix)
    \item C is m x q
    \item D is n x p
    \item c(t) is q x 1
    \item d(t) is p x 1
  }
If a parameter is time-varying then the time dimension is the 3rd dimension. Thus a time-varying Z parameters would be n x m x T where T is the length of the data time series.
  
The main fitting function is \code{\link{MARSS}()} which is used to fit a specified model to data and estimate the model parameters.  \code{\link{MARSS}()} estimates the model parameters using an EM algorithm (primarily but see \code{\link{MARSSoptim}}).  Functions are provided for a variety of output including parameter confidence intervals and the observed Fisher Information matrix, smoothed state estimates with confidence intervals, all the Kalman filter and smoother outputs, residuals and residual diagnostics, printing and plotting, and summaries in tibble format,
}
\details{

\strong{Main MARSS functions}:
\describe{
      \item{\code{\link{MARSS}}}{ Top-level function for specifying and fitting MARSS models. } 
      \item{\code{\link{fitted.marssMLE}}}{ Fitted Y and X estimates as a tibble. Output can be conditioned on all the data (T), data up to t-1, or data up to t.} 
      \item{\code{\link{coef.marssMLE}}}{ Returns the estimated parameters in a variety of formats. }    
      \item{\code{\link{tidy.marssMLE}}}{ Output parameter estimates and fitted Y and X estimates as a tibble. With confidence intervals if requested} 
      \item{\code{\link{augment.marssMLE}}}{ Model residuals as a tibble. } 
      \item{\code{\link{plot.marssMLE}} }{ A series of plots of fits and residuals diagnostics. } 
      \item{\code{\link{autoplot.marssMLE}} }{ A series of plots using ggplot2 of fits and residuals diagnostics. } 
      \item{\code{\link{glance.marssMLE}}}{ Brief summary of fit. } 
      \item{\code{\link{logLik.marssMLE}}}{ Log-likelihood. }
      \item{\code{\link{print.marssMLE}}}{ Prints a wide variety of output from a \code{\link{marssMLE}} oject. }
      \item{\code{\link{print.marssMODEL}}}{ Prints description of the MARSS model. }
      \item{\code{\link{toLatex.marssMODEL}}}{ Outputs a LaTeX version of the model. }
    }

\strong{Other outputs for a fitted model}:
\describe{
      \item{\code{\link{MARSSsimulate}}}{ Produces simulated data from a MARSS model. } 
      \item{\code{\link{MARSSkf}}, \code{\link{MARSSkfas}}, \code{\link{MARSSkfss}} }{ Kalman filters and smoothers with extensive output of all the intermediate filter and smoother variances and expectations. } 
      \item{\code{\link{MARSSaic}}}{ Calculates AICc, AICc, and various bootstrap AICs. } 
      \item{\code{\link{MARSSparamCIs}}}{ Adds confidence intervals to a \code{\link{marssMLE}} object. }
      \item{\code{\link{MARSShessian}}}{ Computes an estimate of the variance-covariance matrix for the MLE parameters. }    
      \item{\code{\link{MARSSFisherI}}}{ Returns the observed Fisher Information matrix. }    
    }

\strong{Important internal MARSS functions (called by the above functions)}:
\describe{
      \item{\code{\link{MARSSkem}}}{ Estimates MARSS parameters using an EM algorithm. } 
      \item{\code{\link{MARSSoptim}}}{ Estimates MARSS parameters using a quasi-Newton algorithm via \code{\link{optim}}. } 
      \item{\code{\link{MARSShatyt}}}{ Calculates the expectations involving Y. } 
      \item{\code{\link{MARSSinnovationsboot}}}{ Creates innovations bootstrapped data. }
      \item{\code{\link{MARSS.marss}}}{ Discusses the form in which MARSS models are stored internally. }
    }
Use \code{help.search("internal", package="MARSS")} to see the documentation of all the internal functions in the MARSS R package.

}
\author{
  Eli Holmes, Eric Ward and Kellie Wills, NOAA, Seattle, USA.

  eli(dot)holmes(at)noaa(dot)gov, eric(dot)ward(at)noaa(dot)gov, 
  
  kellie(dot)wills(at)noaa(dot)gov
}
\references{
The MARSS User Guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   \href{https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf}{User Guide} or type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.

The MARSS Quick Start Guide:  \href{https://cran.r-project.org/package=MARSS/vignettes/Quick_Start.pdf}{Quick Start Guide} or type \code{RShowDoc("Quick_Start",package="MARSS")} to open a copy.

}

\keyword{ package }
