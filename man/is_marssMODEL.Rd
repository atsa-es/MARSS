\name{is.marssMODEL}
\alias{is.marssMODEL}
\alias{is.marssMODEL_marss}
\alias{is.marssMODEL_marxss}
\alias{is.marssMODEL_dfa}
\keyword{internal}

\title{ Test Model Objects }
\description{
These are model objects and utility functions for model objects in the package \code{\link{MARSS-package}}.  
  \code{is.marssMODEL()} ensures model consistency.
  \code{MARSS_formname()} translates a model list as passed in call to \code{\link{MARSS}()} into a marssMODEL model object.
}
\usage{
is.marssMODEL(MODELobj, method="kem")
}

\arguments{
  \item{MODELobj}{ An object of class marssMODEL. }
  \item{method}{ Method used for fitting in case there are special constraints for that method. }
}
\details{
A \code{\link{marssMODEL}} object is an R representation of a MARSS model along with the data.  
    Data in a \code{\link{marssMODEL}} object consists of multivariate time series data in which time is across columns and the n observed time series are in the n different rows.

The base MARSS model (form=marss) is
\describe{
\item{x(t) = B(t) x(t-1) + U(t) + w(t), where w(t) ~ MVN(0,Q(t))}{}
\item{y(t) = Z(t) x(t) + A(t) + v(t), where v(t) ~ MVN(0,R(t))}{}
\item{x(0) ~ MVN(x0, V0) or x(1) ~ MVN(x0, V0)}{} 
}

The marssMODEL(form=marss) object describes this MARSS model but written in vec form:
\describe{
\item{x(t) = kron(x(t-1),I)(f_b(t)+D_b(t)b) + (f_u(t)+D_u(t)u) + w(t), where w(t) ~ MVN(0,Q)}{}
\item{vec(Q) = f_q(t)+D_q(t)q}{} 
\item{y(t) = kron(x(t),I)(f_z(t)+D_z(t)z) + (f_a(t)+D_a(t)a) + v(t), where v(t) ~ MVN(0,R)}{}
\item{vec(R) = f_r(t)+D_r(t)r}{} 
\item{x(0) ~ MVN(f_p+D_p p, V0) or x(1) ~ MVN(f_p+D_p p, V0}{}
\item{vec(V0) = f_l+D_l l}{}  
}
In the marssMODEL(form=marss) object, f(t) + D(t)m, is the vec of a matrix M(t), so f_b(t)+D_b(t)b would be vec(B(t)). The estimated parameters are in the column vectors: b, u, q, z, a, r, p, and l.  Each matrix M(t) is f(t)+D(t)m so is the sum of a fixed part f(t) and the linear combination, D(t), of the free (or estimated) part m.  

The vec form of the MARSS model is specified by 3D matrices for each f and D for each parameter: B, U, Q, Z, A, R, x0, V0.  The number of columns in the D matrix for a parameter determines the number of estimated values for that parameter.  
The first dimension for f (\code{fixed}) and D (\code{free}) must be:
  \describe{
    \item{Z}{ n x m }
    \item{B, Q, and V0}{ m x m }
    \item{U and x0}{ m x 1 }
    \item{A}{ n x 1 }
    \item{R}{ n x n }
  }
  The third dimension of f (\code{fixed}) and D (\code{free}) is either 1 (if not time-varying) or T (if time-varying).   The second dimension of f (\code{fixed}) is always 1, while the second dimension of D  (\code{free}) depends on how many values are being estimated for a matrix.  It can be 0 (if the matrix is fixed) or up to the size of the matrix (if all elements are being estimated).

}
\value{
A vector of error messages or NULL is no errors.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}

\seealso{ \code{\link{MARSS}}, \code{\link{MARSS.marxss}}, \code{\link{marssMODEL}} } 


