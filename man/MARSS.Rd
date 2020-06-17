\name{MARSS}
\alias{MARSS}

\title{ Fit a MARSS Model via Maximum-Likelihood Estimation }
\description{
  This is the main function for fitting multivariate autoregressive state-space (MARSS) models with linear constraints.  Scroll down to the bottom to see some short examples.  To open a guide to show you how to get started quickly, type \code{RShowDoc("Quick_Start",package="MARSS")}.  To open the MARSS User Guide from the command line, type \code{RShowDoc("UserGuide",package="MARSS")}. To get an overview of the package and all its main functions and how to get output (parameter estimates, fitted values, residuals, Kalmin filter or smoother output, or plots), go to \code{\link{MARSS-package}}. If \code{MARSS()} is throwing errors or warnings that you don't understand, try the Troubleshooting section of the user guide or type \code{\link{MARSSinfo}()} at the command line.
  
The default MARSS model form is "marxss", which is Multivariate Auto-Regressive(1) eXogenous inputs State-Space model:
\deqn{\mathbf{x}_{t+1} = \mathbf{B} \mathbf{x}_t + \mathbf{u} + \mathbf{C} \mathbf{c}_t + \mathbf{G} \mathbf{w}_t, \textrm{ where } \mathbf{w}_t \sim \textrm{MVN}(0,\mathbf{Q})}{x(t+1) = B x(t) + u + C c(t) + G w(t), where w(t) ~ MVN(0,Q)}
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{D} \mathbf{d}_t + \mathbf{H} \mathbf{v}_t, \textrm{ where } \mathbf{v}_t \sim \textrm{MVN}(0,\mathbf{R})}{y(t) = Z x(t) + a + D d(t) + H v(t), where v(t) ~ MVN(0,R)}
\deqn{\mathbf{X}_1 \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) \textrm{ or } \mathbf{X}_0 \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) }{X(1) ~ MVN(x0, V0) or X(0) ~  MVN(x0, V0) }
All parameters (except x0 and V0) can be time-varying; \eqn{t} is left off the parameters to remove clutter. All parameters can be zero-ed out, including the variance matrices. Parameters can be estimated or fixed. The parameters are everything except \eqn{mathbf{x}}{x}, \eqn{mathbf{y}}{y}, \eqn{mathbf{v}}{v}, \eqn{mathbf{w}}{w}, \eqn{mathbf{c}}{c} and \eqn{mathbf{c}}{c}. \eqn{mathbf{y}}{y} are data (missing values allowed). c and d are inputs.

The parameters can have fixed values and linear constraints. This is an example of a 3x3 matrix with linear constraints. All matrix elements can be written as a linear function of \eqn{a}, \eqn{b}, and \eqn{c}:
\deqn{\left[\begin{array}{c c c} a+2b & 1 & a\\ 1+3a+b & 0 & b \\ 0 & -2 & c\end{array}\right]}{ [a+2b  1   a \\n 1+3a+b   0   b \\n 0   -2    c ]}
Values such as \eqn{a b} or \eqn{a^2} or \eqn{log(a)} are not linear constraints.

}
\usage{
MARSS(y, 
    model=NULL,
    inits=NULL,
    miss.value=as.numeric(NA),
    method = "kem",
    form = "marxss",
    fit=TRUE, 
    silent = FALSE,
    control = NULL,
    fun.kf = "MARSSkfas",
    ...) 
}

\arguments{
The default settings for the optional arguments are set in \code{MARSSsettings.R} and are given below in the details section. For form specific defaults see the form help file (e.g. \code{\link{MARSS.marxss}} or \code{\link{MARSS.dfa}}).
  \item{y}{ A n x T matrix of n time series over T time steps. Only y is required for the function.}
  \item{inits}{ A list with the same form as the list outputted by \code{coef(fit)} that specifies initial values for the parameters.  See also \code{\link{MARSS.marxss}}. }
  \item{model}{ Model specification using a list of parameter matrix text shortcuts or matrices. See Details and \code{\link{MARSS.marxss}} for the default form. Or better yet open the Quick Start Guide \code{RShowDoc("Quick_Start",package="MARSS")}. }
  \item{miss.value}{ Deprecated.  Denote missing values by NAs in your data. }
  \item{method}{ Estimation method. MARSS provides an EM algorithm (\code{method="kem"}) (see \code{\link{MARSSkem}}) and the BFGS algorithm (\code{method="BFGS"}) (see \code{\link{MARSSoptim}}). }
  \item{form}{ The equation form used in the \code{MARSS()} call.  The default is "marxss". See \code{\link{MARSS.marxss}} or \code{\link{MARSS.dfa}}. }
  \item{fit}{ TRUE/FALSE Whether to fit the model to the data. If FALSE, a \code{\link{marssMLE}} object with only the model is returned. }
  \item{silent}{ TRUE/FALSE Suppresses printing of full error messages, warnings, progress bars and convergence information. Setting silent=2 will produce more verbose error messages and progress information.}
  \item{control}{ Estimation options for the maximization algorithm. The typically used control options for method="kem" are below but see  \code{\link{marssMLE}} for the full list of control options.  Note many of these are not allowed if method="BFGS"; see  \code{\link{MARSSoptim}} for the allowed control options for this method.
    \itemize{
      \item{\code{minit}}{   The minimum number of iterations to do in the maximization routine (if needed by method). If \code{method="kem"}, this is an easy way to up the iterations and see how your estimates are converging. (positive integer)} 
      \item{\code{maxit}}{   Maximum number of iterations to be used in the maximization routine (if needed by method) (positive integer). }
      \item{\code{min.iter.conv.test}}{ Minimum iterations to run before testing convergence via the slope of the log parameter versus log iterations.} 
      \item{\code{conv.test.deltaT=9}}{ Number of iterations to use for the testing convergence via the slope of the log parameter versus log iterations. }
      \item{\code{conv.test.slope.tol}}{ The slope of the log parameter versus log iteration to use as the cut-off for convergence.  The default is 0.5 which is a bit high.  For final analyses, this should be set lower.  If you want to only use abstol as your convergence test, then to something very large, for example conv.test.slope.tol=1000.  Type \code{MARSSinfo(11)} to see some comments of when you might want to do this.}
      \item{\code{abstol}}{ The logLik.(iter-1)-logLik.(iter) convergence tolerance for the maximization routine. To meet convergence both the abstol and slope tests must be passed.} 
      \item{\code{allow.degen}}{ Whether to try setting Q or R elements to zero if they appear to be going to zero. }
      \item{\code{trace}}{ An integer specifying the level of information recorded and error-checking run during the algorithms.  trace=0, specifies basic error-checking and brief error-messages; trace>0 will print full error messages. In addition if trace>0, the Kalman filter output will be added to the outputted \code{marssMLE} object.  Additional information recorded depends on the method of maximization. For the EM algorithm, a record of each parameter estimate for each EM iteration will be added.   See \code{\link{optim}} for trace output details for the BFGS method.   trace=-1 will turn off most internal error-checking and most error messages.  The internal error checks are time expensive so this can speed up MARSS.  This is particularly useful for bootstrapping and simulation studies.}
      \item{\code{silent}}{   TRUE/FALSE(default), Suppresses all printing including progress bars, error messages and convergence information. 0, Turns on all printing of progress bars, fitting information and error messages.  2, Prints a brief success/failure message.    }
      \item{\code{safe}}{   TRUE/FALSE(default), Setting safe=TRUE runs the Kalman smoother after each parameter update rather than running the smoother only once after updated all parameters.  The latter is faster but is not a strictly correct EM algorithm.  In most cases, safe=FALSE (default) will not change the fits.  If this setting does cause problems, you will know because you will see an error regarding the log-likelihood dropping and it will direct you to set safe=TRUE.    }
  } }
    \item{fun.kf}{   What Kalman filter function to use.  MARSS has two: \code{\link{MARSSkfas}()} which is based on the Kalman filter in the \link[KFAS]{KFAS} package based on Koopman and Durbin and \code{\link{MARSSkfss}()} which is a native R implimentation of the Kalman filter and smoother in Shumway and Stoffer. The KFAS filter is much faster.  \code{\link{MARSSkfas}()} modifies the input and output in order to output the lag-one covariance smoother needed for the EM algorithm (per page 321 in Shumway and Stoffer (2000).}
    \item{...}{ Optional arguments passed to function specified by form. }
}
\details{
The \code{model} argument specifies the structure of your model. There is a one-to-one correspondence between how you would write your model in matrix form on the whiteboard and how you specify the model for \code{MARSS()}.  Many different types of multivariate time-series models can be converted to the MARSS form. See the \href{https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf}{User Guide} and \href{https://cran.r-project.org/package=MARSS/vignettes/Quick_Start.pdf}{Quick Start Guide} for examples. 

The MARSS package has two forms for standard users: marxss and dfa.
\describe{
\item{\code{\link{MARSS.marxss}}}{This is the default form. This is a MARSS model with (optional) inputs c(t) or d(t). Most users will want this help page.}
\item{\code{\link{MARSS.dfa}}}{This is a model form to allow easier specification of models for Dynamic Factor Analysis. The Z parameters has a specific form and the Q is set at i.i.d with variance of 1.}
}
Those looking to modify or understand the base code, should look at \code{\link{MARSS.marss}} and 
\code{\link{MARSS.vectorized}}. These describe the forms used by the base functions. The EM algorithm uses the MARSS model written in vectorized form. This form is what allows linear constraints. 

The likelihood surface for MARSS models can be  multimodal or with strong ridges.  It is  recommended that for final analyses the ML estimates are checked by using a Monte Carlo initial conditions search; see the chapter on initial conditions searches in the User Guide.  This requires more computation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.  Also it is wise to check the EM results against the BFGS results (if possible) if there are strong ridges in the likelihood.  Such ridges seems to slow down the EM algorithm considerably and can cause the algorithm to report convergence far from the ML values.  EM steps up the likelihood and the convergence test is based on the rate of change of the LL in each step; once on a strong ridge, the steps can slow dramatically.  You can force the algorithm to keep working by setting \code{minit}.  BFGS seems less hindered by the ridges but can be prodigiously slow for some multivariate problems.  BFGS tends to work better if you give it good initial conditions (see Examples below for how to do this).
 
If you are working with models with time-varying parameters, it is important to notice the time-index for the parameters in the process equation (the x equation).  In some formulations (e.g. in \code{\link[KFAS]{KFAS}}), the process equation is \eqn{x(t)=B(t-1)x(t-1)+w(t-1)} so \eqn{B(t-1)} goes with \eqn{x(t)} not \eqn{B(t)}. Thus one needs to be careful to line up the time indices when passing in time-varying parameters to \code{MARSS()}.  See the User Guide for examples.
}

\value{
  An object of class \code{\link{marssMLE}}.  The structure of this object is discussed below, but if you want to know how to get specific output (like residuals, coefficients, smoothed states, confidence intervals, etc), see \code{\link{print.marssMLE}}, \code{\link{tidy.marssMLE}}, \code{\link{MARSSresiduals}} and \code{\link{plot.marssMLE}}.  
  
  The outputted \code{\link{marssMLE}} object has the following components:
  \item{model}{ MARSS model specification. It is a \code{\link{marssMODEL}} object in the form specified by the user in the \code{MARSS()} call.  This is used by print functions so that the user sees the expected form.}
  \item{marss}{ The \code{\link{marssMODEL}} object in marss form. This form is needed for all the internal algorithms, thus is a required part of a \code{\link{marssMLE}} object.}
  \item{call}{ All the information passed in in the \code{MARSS()} call. }
  \item{start}{ List with specifying initial values that were used for each parameter matrix. } 
  \item{control}{ A list of estimation options, as specified by arguments \code{control}. }
  \item{method}{ Estimation method. }
  If \code{fit=TRUE}, the following are also added to the \code{\link{marssMLE}} object.
  If \code{fit=FALSE}, a \code{\link{marssMLE}} object ready for fitting via the specified \code{method} is returned.
  \item{par}{ A list of estimated parameter values Z, A, R, B, U, Q, x0, V0 in marss form. Use \code{\link{print.marssMLE}}, \code{\link{tidy.marssMLE}} or \code{\link{coef.marssMLE}} for outputing the model estimates in the MARSS() call (e.g. the default "marxss" form).}
  \item{states}{ The expected value of x conditioned on all the data, i.e. smoothed states. }
  \item{states.se}{ The standard errors of the expected value of x. }
  \item{ytT}{ The expected value of y conditioned on the data.  Note this is just y for those y that are not missing. }
  \item{ytT.se}{ The standard errors of the expected value of y. Note this is 0 for any non-missing y.}
  \item{numIter}{ Number of iterations required for convergence. }
  \item{convergence}{ Convergence status. 0 means converged successfully.  Anything else is a warning or error. 2 means the MLEobj has an error; the MLEobj is returned so you can debug it.  The other numbers are errors during fitting.  The error code depends on the fitting method.  See \code{\link{MARSSkem}} and \code{\link{MARSSoptim}}. }
  \item{logLik}{ Log-likelihood. }
  \item{AIC}{ Akaike's Information Criterion. }
  \item{AICc}{ Sample size corrected AIC. }
If \code{control$trace} is set to 1 or greater, the following are also added to the \code{\link{marssMLE}} object.
  \item{kf}{ A list containing Kalman filter/smoother output from \code{\link{MARSSkf}}.  This is not normally added to a \code{\link{marssMLE}} object since it is verbose, but can be added using \code{\link{MARSSkf}()}.}
  \item{Ey}{ A list containing output from \code{\link{MARSShatyt}}. This isn't normally added to a \code{\link{marssMLE}} object since it is verbose, but can be computed using \code{\link{MARSShatyt}()}.}
}
\references{ 
The MARSS User Guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.          

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]

Holmes, E. E., E. J. Ward and K. Wills. (2012) MARSS: Multivariate autoregressive state-space models for analyzing time-series data. R Journal 4: 11-19.
}
\author{ 
  Eli Holmes, Eric Ward and Kellie Wills, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov and eric(dot)ward(at)noaa(dot)gov
}
\seealso{
\code{\link{marssMLE}} \code{\link{MARSSkem}}, \code{\link{MARSSoptim}}, \code{\link{MARSSkf}}, \code{\link{MARSS-package}}, \code{\link{print.marssMLE}}, \code{\link{plot.marssMLE}}, \code{\link{print.marssMODEL}}, \code{\link{MARSS.marxss}}, \code{\link{MARSS.dfa}},  \code{\link{fitted.marssMLE}}, \code{\link{MARSSresiduals}}, \code{\link{predict.marssMLE}},
\code{\link{residuals.marssMLE}}, \code{\link{tidy.marssMLE}}, \code{\link{coef.marssMLE}} 
}
\examples{
dat <- t(harborSealWA)
dat <- dat[2:4, ] # remove the year row
# fit a model with 1 hidden state and 3 observation time series
kemfit <- MARSS(dat, model = list(
  Z = matrix(1, 3, 1),
  R = "diagonal and equal"
))
kemfit$model # This gives a description of the model
print(kemfit$model) # same as kemfit$model
summary(kemfit$model) # This shows the model structure

# add CIs to a marssMLE object
# default uses an estimated Hessian matrix
kem.with.hess.CIs <- MARSSparamCIs(kemfit)
kem.with.hess.CIs

# fit a model with 3 hidden states (default)
kemfit <- MARSS(dat, silent = TRUE) # suppress printing
kemfit

# Fit the above model with BFGS using a short EM fit as initial conditions
kemfit <- MARSS(dat, control=list(minit=5, maxit=5))
bffit <- MARSS(dat, method="BFGS", inits=kemfit)

# fit a model with 3 correlated hidden states
# with one variance and one  covariance
# maxit set low to speed up example, but more iters are needed for convergence
kemfit <- MARSS(dat, model = list(Q = "equalvarcov"), control = list(maxit = 50))
# use Q="unconstrained" to allow different variances and covariances

# fit a model with 3 independent hidden states
# where each observation time series is independent
# the hidden trajectories 2-3 share their U parameter
kemfit <- MARSS(dat, model = list(U = matrix(c("N", "S", "S"), 3, 1)))

# same model, but with fixed independent observation errors
# and the 3rd x processes are forced to have a U=0
# Notice how a list matrix is used to combine fixed and estimated elements
# all parameters can be specified in this way using list matrices
kemfit <- MARSS(dat, model = list(U = matrix(list("N", "N", 0), 3, 1), R = diag(0.01, 3)))

# fit a model with 2 hidden states (north and south)
# where observation time series 1-2 are north and 3 is south
# Make the hidden state process independent with same process var
# Make the observation errors different but independent
# Make the growth parameters (U) the same
# Create a Z matrix as a design matrix that assigns the "N" state to the first 2 rows of dat
# and the "S" state to the 3rd row of data
Z <- matrix(c(1, 1, 0, 0, 0, 1), 3, 2)
# You can use factor is a shortcut making the above design matrix for Z
# Z <- factor(c("N","N","S"))
# name the state vectors
colnames(Z) <- c("N", "S")
kemfit <- MARSS(dat, model = list(
  Z = Z,
  Q = "diagonal and equal", R = "diagonal and unequal", U = "equal"
))

# print the model followed by the marssMLE object
kemfit$model

\dontrun{
# simulate some new data from our fitted model
sim.data <- MARSSsimulate(kemfit, nsim = 10, tSteps = 10)

# Compute bootstrap AIC for the model; this takes a long, long time
kemfit.with.AICb <- MARSSaic(kemfit, output = "AICbp")
kemfit.with.AICb
}

\dontrun{
# Many more short examples can be found in the
# Quick Examples chapter in the User Guide
RShowDoc("UserGuide", package = "MARSS")

# You can find the R scripts from the chapters by
# going to the index page
RShowDoc("index", package = "MARSS")
}

}

