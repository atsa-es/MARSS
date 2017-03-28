\name{MARSS.marxss}
\alias{MARSS.marxss}
\title{ Multivariate AR-1 State-space Model with Inputs }
\description{
The argument \code{form="marxss"} in a MARSS() function call specifies a MAR-1 model with eXogenous variables model.  This is a MARSS(1) model of the form:
\describe{
\item{x(t) = B(t) x(t-1) + u(t) + C(t)c(t) + w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = Z(t) x(t) + a(t) + D(t)d(t) + v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(0) ~ MVN(x0, V0)}{} 
}
Note, marxss is a model form.  A model form is defined by a collection of form functions discussed in \code{\link{marssMODEL}}.  These functions are not exported to the user, but are called by MARSS() using the argument \code{form}.  
}
\section{Usage}{
\code{MARSS(y, 
    inits=NULL,
    model=NULL,
    miss.value=as.numeric(NA),
    method = "kem",
    form = "marxss",
    fit=TRUE, 
    silent = FALSE,
    control = NULL,
    MCbounds = NULL,
    fun.kf = "MARSSkfas",
    ...)}
}

\details{
The allowed arguments when \code{form="marxss"} are 1) the arguments common to all forms: "data", "inits", "MCbounds", "control", "method", "form", "fit", "silent", "fun.kf" (see \code{\link{MARSS}} for information on these arguments) and 2) the argument "model" which is a list describing the MARXSS model (the model list is described below).
See the Quick Start guide (\code{RShowDoc("Quick_Start",package="MARSS")}) or the User Guide (\code{RShowDoc("UserGuide",package="MARSS")}) for examples.

The argument \code{model} must be a list.  The elements in the list specify the structure for the B, u, C, c, Q, Z, a, D, d, R, x0, and V0 in the MARXSS model (above).  The list elements can have the following values:
\itemize{
    \item{\code{Z}} {Default="identity".  A text string, "identity","unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", or "onestate", or a length n vector of factors specifying which of the m hidden state time series correspond to which of the n observation time series.   May be specified as a n x m list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric n x m matrix to use a custom fixed Z. "onestate" gives a n x 1 matrix of 1s. "identity","unconstrained", "diagonal and unequal", "diagonal and equal", and "equalvarcov" all specify n x n matrices.}
    \item{\code{B}} {Default="identity".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". Can also be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use custom fixed B, but in this case all the eigenvalues of B must fall in the unit circle.}
    \item{\code{U}, \code{x0}} {Default="unconstrained". A text string, "unconstrained", "equal", "unequal" or "zero". May be specified as a m x 1 list matrix for general specification of both fixed and shared elements within the matrix. May also be specified as a numeric m x 1 matrix to use a custom fixed U or x0. Notice that U is capitalized.}
     \item{\code{A}} {Default="scaling".  A text string, "scaling","unconstrained", "equal", "unequal" or "zero". May be specified as a  n x 1 list matrix for general specification of both fixed and shared elements within the matrix. May also be specified as a numeric n x 1 matrix to use a custom fixed A. Care must be taken when specifying A so that the model is not under-constrained and unsolveable model.  The default "scaling" only applies to Z matrices that are design matrices (only 1s and 0s and all rows sum to 1). When a column in Z has multiple 1s, the first row with a 1 is assigned A=0 and the rows with 1s for that column have an estimated A. This is used to treat A as an intercept where one A for each X (hidden state) is fixed at 0 and any other Ys associated with that X have an estimated A value.  This ensures a solvable model (when Z is a design matrix). A is capitalized.}
   \item{\code{Q}} {Default="diagonal and unequal".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". May be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use a custom fixed matrix.}
   \item{\code{R}} {Default="diagonal and equal".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". May be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric n x n matrix to use a custom fixed matrix.}
   \item{\code{V0}} {Default="zero".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". May be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric m x m matrix to use a custom fixed matrix.}
    \item{\code{D} and \code{C}} {Default="zero".  A text string, "identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "zero". Can be specified as a list matrix for general specification of both fixed and shared elements within the matrix.  May also be specified as a numeric matrix to use custom fixed values. Must have n rows (D) or m rows (C).}
   \item{\code{d} and \code{c}} {Default="zero".  Numeric matrix.  No missing values allowed.  Must have 1 column or the same number of columns as the data, y.  The numbers of rows in d must be the same as number of columns in D; similarly for c and C.  c and d are lower case.}
   \item{\code{tinitx}} {Default=0.  Whether the initial state is specified at t=0 (default) or t=1.}
  }       
All parameters except x0 and V0 may be time-varying.  If time-varying, then text shortcuts cannot be used.  Enter as an array with the 3rd dimension being time.   Time dimension must be 1 or equal to the number of time-steps in the data.    See Quick Start guide (\code{RShowDoc("Quick_Start",package="MARSS")}) or the User Guide (\code{RShowDoc("UserGuide",package="MARSS")}) for  examples.Valid model structures for \code{method="BFGS"} are the same as for \code{method="kem"}.  See  \code{\link{MARSSoptim}} for the allowed options for this method.

 The default estimation method, \code{method="kem"}, is the EM algorithm described in the user guide.  The default settings for the control and inits arguments are set via \code{MARSS:::alldefaults$kem} in \code{MARSSsettings.R}.  The defaults for the model argument are set in \code{MARSS_marxss.R} For this method, they are:
\itemize{
    \item{inits = list(Z=1, B=1, U=0, Q=0.05, A=0, R=0.05, x0=-99, V0=0)}
    \item{model = list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", Q="diagonal and unequal", x0="unconstrained", V0="zero", C="zero",D="zero",c=matrix(0,0,1), d=matrix(0,0,1), tinitx=0, diffuse=FALSE)}
    \item{control=list(minit=15, maxit=500, abstol=NULL, trace=0,
safe=FALSE, allow.degen=TRUE, min.degen.iter=50, degen.lim=1.0e-04, MCInit=FALSE, numInits = 500, numInitSteps = 10, min.iter.conv.test=15, conv.test.deltaT=9, conv.test.slope.tol= 0.5,
boundsInits=list(B=c(0,1), U=c(-1,1), Q = c(sqrt(0.1),0.1,0.1),
Z=c(0,1), A=c(-1,1), R = c(sqrt(0.1),0.1,0.1) ) )}
  }
For \code{method="BFGS"}, type  \code{MARSS:::alldefaults$BFGS} to see the defaults.


}
\value{
A object of class \code{\link{marssMLE}}.  See \code{\link{print.marssMLE}} for a discussion of the various output available for \code{\link{marssMLE}} objects (coefficients, residuals, Kalman filter and smoother output, imputed values for missing data, etc.).  See \code{\link{MARSSsimulate}} for simulating from \code{\link{marssMLE}} objects. \code{\link{MARSSboot}} for bootstrapping, \code{\link{MARSSaic}} for calculation of various AIC related model selection metrics, and \code{\link{MARSSparamCIs}} for calculation of confidence intervals and bias. 
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  }
\seealso{ \code{\link{marssMODEL}} \code{\link{MARSS.dfa}} } 
\examples{
\dontrun{
#See the MARSS man page for examples
?MARSS

#and the Quick Examples chapter in the User Guide
RShowDoc("UserGuide",package="MARSS")
}
}
