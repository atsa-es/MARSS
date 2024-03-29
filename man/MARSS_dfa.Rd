\name{MARSS.dfa}
\alias{MARSS.dfa}

\title{ Multivariate Dynamic Factor Analysis }

\description{
The Dynamic Factor Analysis model in MARSS is
The argument \code{form="marxss"} in a \code{\link{MARSS}()} function call specifies a MAR-1 model with eXogenous variables model.  This is a MARSS(1) model of the form:
\deqn{\mathbf{x}_{t} = \mathbf{x}_{t-1} + \mathbf{w}_t, \textrm{ where } \mathbf{W}_t \sim \textrm{MVN}(0,\mathbf{I})}{x(t) = x(t-1) + w(t), where W(t) ~ MVN(0,I)}
\deqn{\mathbf{y}_t = \mathbf{Z}_t \mathbf{x}_t + \mathbf{D}_t \mathbf{d}_t + \mathbf{v}_t, \textrm{ where } \mathbf{V}_t \sim \textrm{MVN}(0,\mathbf{R}_t)}{y(t) = Z(t) x(t) + D(t) d(t) + v(t), where V(t) ~ MVN(0,R(t))}
\deqn{\mathbf{X}_1 \sim \textrm{MVN}(\mathbf{x0}, 5\mathbf{I})}{X(1) ~ MVN(x0, 5I) }
Note, by default \eqn{\mathbf{x}_1}{x(1)} is treated as a diffuse prior.

Passing in \code{form="dfa"} to \code{\link{MARSS}()} invokes a helper function to create that model and creates the \eqn{\mathbf{Z}}{Z} matrix for the user.  \eqn{\mathbf{Q}}{Q} is by definition identity, \eqn{\mathbf{x}_0}{x0} is zero and \eqn{\mathbf{V_0}}{V0} is diagonal with large variance (5). \eqn{\mathbf{u}}{U} is zero, \eqn{\mathbf{a}}{A} is zero, and covariates only enter the \eqn{\mathbf{y}}{Y} equation.  Because \eqn{\mathbf{u}}{U} and \eqn{\mathbf{a}}{A} are 0, the data should have mean 0 (demeaned) otherwise one is likely to be creating a structurally inadequate model (i.e. the model implies that the data have mean = 0, yet data do not have mean = 0 ).
}
\section{Usage}{
\code{MARSS(y, 
    inits = NULL, 
    model = NULL, 
    miss.value = as.numeric(NA), 
    method = "kem", 
    form = "dfa", 
    fit = TRUE, 
    silent = FALSE, 
    control = NULL, 
    fun.kf = "MARSSkfas", 
    demean = TRUE, 
    z.score = TRUE)}
}

\arguments{
Some arguments are common to all forms: "y" (data), "inits", "control", "method", "form", "fit", "silent", "fun.kf". See \code{\link{MARSS}} for information on these arguments.

In addition to these, form="dfa" has some special arguments that can be passed in:
\itemize{
    \item \code{demean}  Logical.  Default is TRUE, which means the data will be demeaned.
    \item \code{z.score}  Logical.  Default is TRUE, which means the data will be z-scored (demeaned and variance standardized to 1).
    \item \code{covariates}  Covariates (\eqn{d}) for the \eqn{y} equation. No missing values allowed and must be a matrix with the same number of time steps as the data. An unconstrained \eqn{D} matrix will estimated.
    }

The \code{model} argument of the \code{\link{MARSS}()} call is constrained in terms of what parameters can be changed and how they can be changed.  See details below. An additional element, \code{m}, can be passed into the \code{model} argument that specifies the number of hidden state variables.  It is not necessarily for the user to specify \code{Z} as the helper function will create a \code{Z} appropriate for a DFA model.
}
\details{
The \code{model} argument is a list. The following details what list elements can be passed in:
    \itemize{
    \item \code{B}  "Identity". The standard (and default) DFA model has B="identity". However it can be "identity", "diagonal and equal", "diagonal and unequal"  or a time-varying fixed or estimated diagonal matrix.
    \item \code{U}  "Zero". Cannot be changed or passed in via model argument.
    \item \code{Q}  "Identity". The standard (and default) DFA model has Q="identity". However, it can be "identity", "diagonal and equal", "diagonal and unequal" or a time-varying fixed or estimated diagonal matrix.
    \item \code{Z}  Can be passed in as a (list) matrix if the user does not want a default DFA \code{Z} matrix.  There are many equivalent ways to construct a DFA \code{Z} matrix.  The default is Zuur et al.'s form (see User Guide).
    \item \code{A}  Default="zero". Can be "unequal", "zero" or a matrix.
    \item \code{R}  Default="diagonal and equal". Can be set to "identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", or a (list) matrix to specify general forms.
    \item \code{x0} Default="zero". Can be "unconstrained", "unequal", "zero", or a (list) matrix.
    \item \code{V0} Default=diagonal matrix with 5 on the diagonal. Can be "identity", "zero", or a matrix.
    \item \code{tinitx} Default=0. Can be 0 or 1.  Tells MARSS whether x0 is at t=0 or t=1.
    \item \code{m} Default=1. Can be 1 to n (the number of y time-series). Must be integer.
    }

See the \href{https://cran.r-project.org/package=MARSS/vignettes/UserGuide.pdf}{User Guide} chapter on Dynamic Factor Analysis for examples of of using \code{form="dfa"}.
}
\value{
A object of class \code{\link{marssMLE}}.  See \code{\link[=print.marssMLE]{print}()} for a discussion of the various output available for \code{\link{marssMLE}} objects (coefficients, residuals, Kalman filter and smoother output, imputed values for missing data, etc.).  See \code{\link{MARSSsimulate}()} for simulating from \code{\link{marssMLE}} objects. \code{\link{MARSSboot}()} for bootstrapping, \code{\link{MARSSaic}()} for calculation of various AIC related model selection metrics, and \code{\link{MARSSparamCIs}()} for calculation of confidence intervals and bias. 
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA. 
  }

\seealso{ \code{\link{MARSS}()}, \code{\link{MARSS.marxss}()} } 

\references{ 
The MARSS User Guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science
Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.          
}
\examples{
\dontrun{
dat <- t(harborSealWA[,-1])
# DFA with 3 states; used BFGS because it fits much faster for this model
fit <- MARSS(dat, model = list(m=3), form="dfa", method="BFGS")

# See the Dynamic Factor Analysis chapter in the User Guide
RShowDoc("UserGuide", package = "MARSS")
}
}

