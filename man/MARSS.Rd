\name{MARSS}
\alias{MARSS}
\alias{MARSSsettings}
\title{ MARSS Model Specification and Estimation }
\description{
  This is the main MARSS function for fitting multivariate autoregressive state-space (MARSS) models.  Scroll down to the bottom to see some short examples.  To open the user guide from the command line, type \code{RShowDoc("UserGuide",package="MARSS")}. To open a guide to show you how to get started quickly, type \code{RShowDoc("Quick_Start",package="MARSS")}.  To open an overview page with package information and how to find all the R code for the user guide chapters, type \code{RShowDoc("index",package="MARSS")}. To get info on the axillary functions (like for bootstrapping and confidence intervals) go to \code{\link{MARSS-package}}.  To see a discussion of how to get output from your model fits, go to \code{\link{print.MARSS}}. If MARSS() is throwing errors or warnings that you don't understand, try the Troubleshooting section of the user guide or type \code{MARSSinfo()} at the command line.
  
The background section on this page is focused on fitting MARSS models in vectorized form. This form will almost certainly look unfamilar.  MARSS works by converting the users' (more familiar model form) into the vectorized form which allows general linear constraints.  You should go to the help page for the form of the model you are fitting to get background on that model form.  Currently the MARSS package has two model forms: marxss and dfa.
\describe{
\item{\code{\link{MARSS.marxss}}}{This is the default form. X_t = B_t*X_t-1 + U_t + C_t*c_t + w_t; Y_t=Z_t*X_t + A_t + D_t*d_t + v_t. Any parameters can be set to zero. Most users will want this help page.}
\item{\code{\link{MARSS.dfa}}}{This is a model form to allow easier specification of models for Dynamic Factor Analysis. X_t = X_t-1 + U_t + w_t; Y_t=Z_t*X_t + A_t + D_t*d_t + v_t. The Z parameters has a specific form and the Q is set at i.i.d with variance of 1.}
}
The rest of this help page discusses the vectorized form of a MARSS model.

The MARSS package fits time-varying state-space models that can be transformed into the form (termed form=marss):
\describe{
\item{x(t) = (t(x(t-1)) (x) I_m)(f_b(t)+D_b(t)beta) + (f_u(t)+D_u(t)upsilon) + w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = (t(x(t)) (x) I_n)(f_z(t)+D_z(t)zeta) + (f_a(t)+D_a(t)alpha) + v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0)}{} 
}
where beta, upsilon, zeta, and alpha are column vectors of estimated values, the f are column vectors of inputs, and the D are matrices of inputs.  The f and D are potentially time-varying. (x) means kronecker product and I_p is a p x p identity matrix.  The function \code{MARSS()} is used to fit MARSS models using the argument \code{form} to specify the type of state-space model being fit.  

Most commonly used multivariate autoregressive state-space models can be reformulated into the form above.  The user is not required to specify their model in the marss form (which is unfamiliar and unwieldy).  Instead \code{MARSS()} uses the \code{form} argument to specify a more familiar state-space form.  The user specifies their model in that (more familiar) form.  \code{MARSS()} calls a helper function \code{MARSS_form} to translate the user's model into form=marss.

The default MARSS form is "marxss" which is the state-space model:
\describe{
\item{x(t) = B(t) x(t-1) + u(t) + C(t)c(t) + G(t)w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = Z(t) x(t) + a(t) + D(t)d(t) + H(t)v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0)}{} 
}
See \code{\link{MARSS.marxss}} for arguments and defaults information.

If you are working with models with time-varying parameters, it is important to notice the time-index for the parameters in the process equation (the x equation).  In some formulations (e.g. in the \code{\link[KFAS]{KFAS}}), the process equation is x(t)=B(t-1)x(t-1)+w(t-1) so B(t-1) goes with x(t) not B(t). Thus one needs to be careful to line up the time indices when passing in time-varying parameters to \code{MARSS()}.
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
  \item{y}{ A n x T matrix of n time series over T time steps. }
  \item{inits}{ A list with the same form as the list outputted by \code{coef(fit)} that specifies initial values for the parameters.  See also \code{\link{MARSS.marxss}}. }
  \item{model}{ Model specification using parameter model text shortcuts or matrices. See Details and \code{\link{MARSS.marxss}} for the default form. Or better yet open the Quick Start Guide \code{RShowDoc("Quick_Start",package="MARSS")}. }
  \item{miss.value}{ Deprecated.  Denote missing values by NAs in your data. }
  \item{method}{ Estimation method. MARSS provides an EM algorithm (\code{method="kem"}) (see \code{\link{MARSSkem}}) and the BFGS algorithm (\code{method="BFGS"}) (see \code{\link{MARSSoptim}}). }
  \item{form}{ The equation form used in the \code{MARSS()} call.  The default is "marxss". See \code{\link{MARSS.marxss}} or \code{\link{MARSS.dfa}}. }
  \item{fit}{ TRUE/FALSE Whether to fit the model to the data. If FALSE, a marssMLE object with only the model is returned. }
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
    \item{fun.kf}{   What Kalman filter function to use.  MARSS has two: MARSSkfas which is based on the Kalman filter in the KFAS package based on Koopman and Durbin and MARSSkfss which is a native R implimentation of the Kalman filter and smoother in Shumway and Stoffer. The KFAS filter is much faster.  MARSSkfas modifies the input and output in order to output the lag-one covariance smoother needed for the EM algorithm (per page 321 in Shumway and Stoffer (2000).}
    \item{...}{ Optional arguments passed to function specified by form. }
}
\details{
 MARSS provides an interface to the base \code{\link{MARSS-package}} functions and allows specification and fitting of MARSS models. The available estimation methods are maximum-likelihood via an EM algorithm (\code{method="kem"}) or via a quasi-Newton algorithm provided by function \code{optim} (\code{method="BFGS"}).  The function \code{MARSS()} allows the user to specify models using the \code{model} argument.  See \code{\link{MARSS.marxss}} for the format of the \code{model} argument for the default \code{marxss} form.  See also the User Guide (reference and link below) or Quick Start Guide.
 
A call to \code{MARSS()} returns an object of class \code{marssMLE}.  The MARSS package has print, coef, and residuals functions that will handle \code{marssMLE} objects.  See \code{\link{print.marssMLE}}, \code{\link{coef.marssMLE}}, and \code{\link{residuals.marssMLE}}. The help page for \code{\link{print.marssMLE}} summarizes all the different output available for \code{marssMLE} objects and describes what the output is mathematically.  Thus, this page is the first place to start and will direct you to the appropriate other method functions (like \code{coef}).

 Many different types of multivariate time-series models can be converted to the MARSS form (see the User Guide). \code{MARSS()} allows the user to specify the form of their model using the argument \code{form}. The default form is a multivariate lag-1 (time-varying) state-space model which is denoted using \code{form="marxss"}. Look at \code{\link{MARSS.marxss}}) to see how the model argument for the marxss form is specified.   The model argument is a list, but the elements of that and how it is converted to a marssMODEL(form=marss) object (needed for the internal algorithms) depends intimately on the equation form.  Thus you will need to refer to that appropriate help page (\code{MARSS.formname}) for your equation form.  

 The likelihood surface for MARSS models can be  multimodal or with strong ridges.  It is  recommended that for final analyses the ML estimates are checked by using a Monte Carlo initial conditions search; see the chapter on initial conditions searches in the User Guide.  This requires more computation time, but reduces the chance of the algorithm terminating at a local maximum and not reaching the true MLEs.  Also it is wise to check the EM results against the BFGS results (if possible) since if there are strong ridges in the likelihood.  Such ridges seems to slow down the EM algorithm considerably and can cause the algorithm to report convergence far from the ML values.  EM steps up the likelihood and the convergence test is based on the rate of change of the LL in each step; once on a strong ridge, the steps can slow dramatically.  You can force the algorithm to keep working by setting \code{minit}.  BFGS seems less hindered by the ridges but can be prodigiously slow for some multivariate problems.
}

\value{
  An object of class \code{\link{marssMLE}}.  The structure of this object is discussed below, but if you want to know how to get specific output (like residuals, coefficients, smoothed states, confidence intervals, etc), go here \code{\link{print.MARSS}}.  
  
  The outputted \code{\link{marssMLE}} object has the following components:
  \item{model}{ MARSS model specification. It is a \code{marssMODEL} object in the form specified by the user in the \code{MARSS()} call.  This is used by print functions so that the user sees the expected form.}
  \item{marss}{ The \code{marssMODEL} object in marss form. This form is needed for all the internal algorithms, thus is a required part of a marssMLE object.}
  \item{call}{ All the information passed in in the \code{MARSS()} call. }
  \item{start}{ List with specifying initial values that were used for each parameter matrix. } 
  \item{control}{ A list of estimation options, as specified by arguments \code{control}. }
  \item{method}{ Estimation method. }
  If \code{fit=TRUE}, the following are also added to the marssMLE object.
  If \code{fit=FALSE}, an marssMLE object ready for fitting via the specified \code{method} is returned.
  \item{par}{ A list of estimated parameter values Z, A, R, B, U, Q, x0, V0. See \code{\link{print.marssMLE}} or \code{\link{coef.marssMLE}} for information on outputing the model estimates.  This will be in form "marss".  Use \code{print} or \code{coef} to output the estimated parameters in the form in the MARSS() call (e.g. the default "marxss" form).}
  \item{states}{ The expected value of x conditioned on the data. }
  \item{states.se}{ The standard errors of the expected value of x. }
  \item{ytT}{ The expected value of y conditioned on the data.  Note this is just y for those y that are not missing. }
  \item{ytT.se}{ The standard errors of the expected value of y. Note this is 0 for any non-missing y.}
  \item{numIter}{ Number of iterations required for convergence. }
  \item{convergence}{ Convergence status. 0 means converged successfully.  Anything else is a warning or error. 2 means the MLEobj has an error; the MLEobj is returned so you can debug it.  The other numbers are errors during fitting.  The error code depends on the fitting method.  See \code{\link{MARSSkem}} and \code{\link{MARSSoptim}}. }
  \item{logLik}{ Log-likelihood. }
  \item{AIC}{ Akaike's Information Criterion. }
  \item{AICc}{ Sample size corrected AIC. }
If \code{control$trace} is set to 1 or greater, the following are also added to the marssMLE object.
  \item{kf}{ A list containing Kalman filter/smoother output from \code{\link{MARSSkf}}.  This isn't normally added to a marssMLE object since it is verbose, but can be computed using \code{MARSSkf(marssMLE)}.}
  \item{Ey}{ A list containing output from \code{\link{MARSShatyt}}. This isn't normally added to a marssMLE object since it is verbose, but can be computed using \code{MARSShatyt(marssMLE)}.}
}
\references{ 
The user guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.          

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]

Holmes, E. E., E. J. Ward and K. Wills. (2012) MARSS: Multivariate autoregressive state-space models for analyzing time-series data. R Journal 4: 11-19.
}
\author{ 
  Eli Holmes, Eric Ward and Kellie Wills, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov, kellie(dot)wills(at)noaa(dot)gov
}
\seealso{
\code{\link{marssMLE}} \code{\link{MARSSkem}} \code{\link{MARSSoptim}} \code{\link{MARSS-package}} \code{\link{print.MARSS}} \code{\link{MARSS.marxss}} \code{\link{MARSS.dfa}}
\code{\link{augment.marssMLE}} \code{\link{tidy.marssMLE}} \code{\link{coef.marssMLE}} 
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
kemfit

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

