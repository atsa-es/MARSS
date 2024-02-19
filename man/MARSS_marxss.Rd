\name{MARSS.marxss}
\alias{MARSS.marxss}
\keyword{appendix}

\title{ Multivariate AR-1 State-space Model with Inputs }
\description{
The argument \code{form="marxss"} in a \code{\link{MARSS}()} function call specifies a MAR-1 model with eXogenous variables model.  This is a MARSS(1) model of the form:
\deqn{\mathbf{x}_{t} = \mathbf{B}_t \mathbf{x}_{t-1} + \mathbf{u}_t + \mathbf{C}_t \mathbf{c}_t + \mathbf{G}_t \mathbf{w}_t, \textrm{ where } \mathbf{W}_t \sim \textrm{MVN}(0,\mathbf{Q}_t)}{x(t) = B(t) x(t-1) + u(t) + C(t) c(t) + G(t) w(t), where W(t) ~ MVN(0,Q(t))}
\deqn{\mathbf{y}_t = \mathbf{Z}_t \mathbf{x}_t + \mathbf{a}_t + \mathbf{D}_t \mathbf{d}_t + \mathbf{H}_t \mathbf{v}_t, \textrm{ where } \mathbf{V}_t \sim \textrm{MVN}(0,\mathbf{R}_t)}{y(t) = Z(t) x(t) + a(t) + D(t) d(t) + H(t) v(t), where V(t) ~ MVN(0,R(t))}
\deqn{\mathbf{X}_1 \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) \textrm{ or } \mathbf{X}_0 \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) }{X(1) ~ MVN(x0, V0) or X(0) ~  MVN(x0, V0) }
Note, by default \eqn{\mathbf{V0}}{V0} is a matrix of all zeros and thus \eqn{\mathbf{x}_1}{x(1)} or \eqn{\mathbf{x}_0}{x(0)} is treated as an estimated parameter not a diffuse prior.

Note, "marxss" is a model form.  A model form is defined by a collection of form functions discussed in \code{\link{marssMODEL}}.  These functions are not exported to the user, but are called by \code{\link{MARSS}()} using the argument \code{form}.  
}
\section{Usage}{
\code{MARSS(y, 
    inits = NULL, 
    model = NULL, 
    miss.value = as.numeric(NA), 
    method = "kem", 
    form = "marxss", 
    fit = TRUE, 
    silent = FALSE, 
    control = NULL, 
    fun.kf = "MARSSkfas", 
    ...)}
}

\details{
The allowed arguments when \code{form="marxss"} are 1) the arguments common to all forms: "y" (data), "inits", "control", "method", "form", "fit", "silent", "fun.kf" (see \code{\link{MARSS}} for information on these arguments) and 2) the argument "model" which is a list describing the MARXSS model (the model list is described below).
See the \href{https://cran.r-project.org/package=MARSS/vignettes/Quick_Start.html}{Quick Start Guide} guide or the \href{https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf}{User Guide} for examples.

The argument \code{model} must be a list.  The elements in the list specify the structure for the \eqn{\mathbf{B}}{B}, \eqn{\mathbf{u}}{u}, \eqn{\mathbf{C}}{C}, \eqn{\mathbf{c}}{c}, \eqn{\mathbf{Q}}{Q}, \eqn{\mathbf{Z}}{Z}, \eqn{\mathbf{a}}{a}, \eqn{\mathbf{D}}{D}, \eqn{\mathbf{d}}{d}, \eqn{\mathbf{R}}{R}, \eqn{\mathbf{x}_0}{x0}, and \eqn{\mathbf{V}_0}{V0} in the MARXSS model (above).  The list elements can have the following values:
\describe{
    \item{\code{Z}}{Default="identity".  A text string, "identity","unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", or "onestate", or a length n vector of factors specifying which of the m hidden state time series correspond to which of the n observation time series.   May be specified as a n x m list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric n x m matrix to use a custom fixed \eqn{\mathbf{Z}}{Z}. "onestate" gives a n x 1 matrix of 1s. "identity","unconstrained", "diagonal and unequal", "diagonal and equal", and "equalvarcov" all specify n x n matrices.}
    \item{\code{B}}{Default="identity".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". Can also be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use custom fixed \eqn{\mathbf{B}}{B}, but in this case all the eigenvalues of \eqn{\mathbf{B}}{B} must fall in the unit circle.}
    \item{\code{U}, \code{x0}}{Default="unconstrained". A text string, "unconstrained", "equal", "unequal" or "zero". May be specified as a m x 1 list matrix for general specification of both fixed and shared elements within the matrix. May also be specified as a numeric m x 1 matrix to use a custom fixed \eqn{\mathbf{u}}{u} or \eqn{\mathbf{x}_0}{x(0)}. Notice that \code{U} is capitalized in the \code{model} argument and output lists.}
     \item{\code{A}}{Default="scaling".  A text string, "scaling","unconstrained", "equal", "unequal" or "zero". May be specified as a  n x 1 list matrix for general specification of both fixed and shared elements within the matrix. May also be specified as a numeric n x 1 matrix to use a custom fixed \eqn{\mathbf{a}}{a}. Care must be taken when specifying \code{A} so that the model is not under-constrained and unsolvable model.  The default, "scaling", only applies to \eqn{\mathbf{Z}}{Z} matrices that are design matrices (only 1s and 0s and all rows sum to 1). When a column in \eqn{\mathbf{Z}}{Z} has multiple 1s, the first row in the \eqn{\mathbf{a}}{a} matrix associated with those \eqn{\mathbf{Z}}{Z} rows is 0 and the other associated \eqn{\mathbf{a}}{a} rows have an estimated value. This is used to treat \eqn{\mathbf{a}}{a} as an intercept where one intercept for each \eqn{\mathbf{x}}{x} (hidden state) is fixed at 0 and any other intercepts associated with that \eqn{\mathbf{x}}{x} have an estimated intercept.  This ensures a solvable model when \eqn{\mathbf{Z}}{Z} is a design matrix. Note in the model argument and output, \code{A} is capitalized.}
   \item{\code{Q}}{Default="diagonal and unequal".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". May be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric g x g matrix to use a custom fixed matrix. Default value of g is m, so \eqn{\mathbf{Q}}{Q} is a m x m matrix.  g is the number of columns in \eqn{\mathbf{G}}{G} (below).}
   \item{\code{R}}{Default="diagonal and equal".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". May be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric h x h matrix to use a custom fixed matrix. Default value of h is n, so \eqn{\mathbf{R}}{R} is a n x n matrix.  h is the num of columns in \eqn{\mathbf{H}}{H} (below).}
   \item{\code{V0}}{Default="zero".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". May be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use a custom fixed matrix.}
    \item{\code{D} and \code{C}}{Default="zero".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". Can be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric matrix to use custom fixed values. Must have n rows (\eqn{\mathbf{D}}{D}) or m rows (\eqn{\mathbf{C}}{C}).}
   \item{\code{d} and \code{c}}{Default="zero".  Numeric matrix.  No missing values allowed.  Must have 1 column or the same number of columns as the data, \eqn{\mathbf{y}}{y}.  The numbers of rows in \eqn{\mathbf{d}}{d} must be the same as number of columns in \eqn{\mathbf{D}}{D}; similarly for \eqn{\mathbf{c}}{c} and \eqn{\mathbf{C}}{C}.}
    \item{\code{G} and \code{H}}{Default="identity".  A text string, "identity". Can be specified as a numeric matrix or array for time-varying cases. Must have m rows and g columns (\eqn{\mathbf{G}}{G}) or n rows and h columns (\eqn{\mathbf{H}}{H}). g is the dim of \eqn{\mathbf{Q}}{Q} and h is the dim of \eqn{\mathbf{R}}{R}.}
   \item{\code{tinitx}}{Default=0.  Whether the initial state is specified at t=0 (default) or t=1.}
  }       
All parameters except \eqn{\mathbf{x}_0}{x0} and \eqn{\mathbf{V}_0}{V0} may be time-varying.  If time-varying, then text shortcuts cannot be used.  Enter as an array with the 3rd dimension being time.   Time dimension must be 1 or equal to the number of time-steps in the data.    See Quick Start guide (\code{RShowDoc("Quick_Start",package="MARSS")}) or the User Guide (\code{RShowDoc("UserGuide",package="MARSS")}) for  examples.Valid model structures for \code{method="BFGS"} are the same as for \code{method="kem"}.  See  \code{\link{MARSSoptim}()} for the allowed options for this method.

 The default estimation method, \code{method="kem"}, is the EM algorithm described in the MARSS User Guide.  The default settings for the control and inits arguments are set via \code{MARSS:::alldefaults$kem} in \code{MARSSsettings.R}.  The defaults for the model argument are set in \code{MARSS_marxss.R} For this method, they are:
\itemize{
    \item{inits = list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=0.05, G=0, H=0, L=0, C=0, D=0, c=0, d=0)}
    \item{model = list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained", V0="zero", C="zero",D="zero",c=matrix(0,0,1), d=matrix(0,0,1), tinitx=0, diffuse=FALSE)}
    \item{control=list(minit=15, maxit=500, abstol=0.001, trace=0, sparse=FALSE,
                   safe=FALSE, allow.degen=TRUE, min.degen.iter=50, degen.lim=1.0e-04, 
                   min.iter.conv.test=15, conv.test.deltaT=9, conv.test.slope.tol= 0.5,  demean.states=FALSE) You can read about these in \code{\link{MARSS}()}.  If you want to speed up your fits, you can turn off most of the model checking using \code{trace=-1}. }
    \item{fun.kf = "MARSSkfas"; This sets the Kalman filter function to use.  \code{MARSSkfas()} is generally more stable as it uses Durban & Koopman's algorithm.  But it may dramatically slow down when the data set is large (more than 10 rows of data).  Try the classic Kalman filter algorithm to see if it runs faster by setting \code{fun.kf="MARSSkfss"}. You can read about the two algorithms in \code{\link{MARSSkf}}.}
  }
For \code{method="BFGS"}, type  \code{MARSS:::alldefaults$BFGS} to see the defaults.


}
\value{
A object of class \code{\link{marssMLE}}.  See \code{\link{print.marssMLE}} for a discussion of the various output available for \code{\link{marssMLE}} objects (coefficients, residuals, Kalman filter and smoother output, imputed values for missing data, etc.).  See \code{\link{MARSSsimulate}} for simulating from \code{\link{marssMLE}} objects. \code{\link{MARSSboot}} for bootstrapping, \code{\link{MARSSaic}} for calculation of various AIC related model selection metrics, and \code{\link{MARSSparamCIs}} for calculation of confidence intervals and bias. See \code{\link{plot.marssMLE}} for some default plots of a model fit.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  }
\seealso{ \code{\link{marssMODEL}}, \code{\link{MARSS.dfa}()} } 
\examples{
\dontrun{
#See the MARSS man page for examples
?MARSS

#and the Quick Examples chapter in the User Guide
RShowDoc("UserGuide",package="MARSS")
}
}
