\name{MARSS.marss}
\alias{MARSS.marss}
\keyword{internal}

\title{ Multivariate AR-1 State-space Model  }
\description{
The form of MARSS models for users is marxss, the MARSS models with inputs.  See \code{\link{MARSS.marxss}}.  In the internal algorithms (e.g. \code{\link{MARSSkem}}), the "marss" form is used and the D and d are incorporated into the A matrix and C and c are incorporated into the U matrix.

This is a MARSS(1) model of the marss form:
\deqn{x(t+1) = B(t) x(t) + U(t)  + G(t) w(t), where w(t) ~ MVN(0,Q(t))}
\deqn{y(t) = Z(t) x(t) + A(t) + H(t) v(t), where v(t) ~ MVN(0,R(t))}
\deqn{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0) }
Note, by default \eqn{V0} is a matrix of all zeros and thus \eqn{x(1)} or \eqn{x(0)} is treated as an estimated parameter not a diffuse prior.


Note, "marss" is a model form.  A model form is defined by a collection of form functions discussed in \code{\link{marssMODEL}}.  These functions are not exported to the user, but are called by \code{\link{MARSS}()} using the argument \code{form}.  These internal functions convert the users model list into the vec form of a MARSS model and do extensive error-checking.
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
    fun.kf = "MARSSkfas",
    ...)}
}

\details{
See the help page for the \code{\link{MARSS.marxss}} form for details. 
}
\value{
A object of class \code{\link{marssMLE}}.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  }
\seealso{ \code{\link{marssMODEL}} \code{\link{MARSS.marxss}} } 
\examples{
\dontrun{
# See the MARSS man page for examples
?MARSS

# and the Quick Examples chapter in the User Guide
RShowDoc("UserGuide", package = "MARSS")
}
}
