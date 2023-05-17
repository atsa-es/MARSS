\name{MARSS.marss}
\alias{MARSS.marss}
\keyword{appendix}


\title{ Multivariate AR-1 State-space Model  }
\description{
The form of MARSS models for users is "marxss", the MARSS models with inputs.  See \code{\link{MARSS.marxss}}.  In the internal algorithms (e.g. \code{\link{MARSSkem}}), the "marss" form is used and the \eqn{\mathbf{D}\mathbf{d}_t}{Dd(t)} are incorporated into the \eqn{\mathbf{a}_t}{a(t)} matrix and \eqn{\mathbf{C}\mathbf{c}_t}{Cc(t)} are incorporated into the \eqn{\mathbf{u}_t}{u(t)}. The \eqn{\mathbf{a}}{a} and \eqn{\mathbf{u}}{u} matrices then become time-varying if the model includes \eqn{\mathbf{d}_t}{d(t)} and \eqn{\mathbf{c}_t}{c(t)}.

This is a MARSS(1) model of the marss form:
\deqn{\mathbf{x}_{t} = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u}_t + \mathbf{G} \mathbf{w}_t, \textrm{ where } \mathbf{w}_t \sim \textrm{MVN}(0,\mathbf{Q})}{x(t) = B x(t-1) + u(t) + G w(t), where w(t) ~ MVN(0,Q)}
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a}_t + \mathbf{H} \mathbf{v}_t, \textrm{ where } \mathbf{v}_t \sim \textrm{MVN}(0,\mathbf{R})}{y(t) = Z x(t) + a(t) + H v(t), where v(t) ~ MVN(0,R)}
\deqn{\mathbf{X}_1 \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) \textrm{ or } \mathbf{X}_0 \sim \textrm{MVN}(\mathbf{x0}, \mathbf{V0}) }{X(1) ~ MVN(x0, V0) or X(0) ~  MVN(x0, V0) }Note, by default \eqn{\mathbf{V0}}{V0} is a matrix of all zeros and thus \eqn{\mathbf{x}_1}{x(1)} or \eqn{\mathbf{x}_0}{x(0)} is treated as an estimated parameter not a diffuse prior. To remove clutter, the rest of the parameters are shown as time-constant (no \eqn{t} subscript) but all parameters can be time-varying.


Note, "marss" is a model form.  A model form is defined by a collection of form functions discussed in \code{\link{marssMODEL}}.  These functions are not exported to the user, but are called by \code{\link{MARSS}()} using the argument \code{form}.  These internal functions convert the users model list into the vec form of a MARSS model and do extensive error-checking.
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
See the help page for the \code{\link{MARSS.marxss}} form for details. 
}
\value{
A object of class \code{\link{marssMLE}}.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.

  }
\seealso{ \code{\link{marssMODEL}}, \code{\link{MARSS.marxss}()} } 
\examples{
\dontrun{
# See the MARSS man page for examples
?MARSS

# and the Quick Examples chapter in the User Guide
RShowDoc("UserGuide", package = "MARSS")
}
}
