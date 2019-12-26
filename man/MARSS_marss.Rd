\name{MARSS.marss}
\alias{MARSS.marss}
\title{ Multivariate AR-1 State-space Model  }
\description{
The form of MARSS models for users is marxss, the MARSS models with inputs.  See \code{\link{MARSS.marxss}}.  In the internal algorithms, the marss form is used and the D and d are incorporated into the A matrix and C and c are incorporated into the U matrix.

This is a MARSS(1) model of the marss form:
\describe{
\item{x(t) = B(t) x(t-1) + u(t) + G(t)w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = Z(t) x(t) + a(t) + H(t)v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(0) ~ MVN(x0, V0)}{} 
}

Note, marss is a model form.  A model form is defined by a collection of form functions discussed in \code{\link{marssMODEL}}.  These functions are not exported to the user, but are called by MARSS() using the argument \code{form}.  These internal functions convert the users model list into the vec form of a MARSS model and do extensive error-checking.
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
