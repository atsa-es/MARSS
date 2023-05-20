\name{MARSS.vectorized}
\alias{MARSS.vectorized}
\keyword{experimental}

\title{ Vectorized Multivariate AR-1 State-space Model  }

\description{
The EM algorithm (\code{\link{MARSSkem}}) in the MARSS package works by converting the more familiar MARSS model in matrix form into the vectorized form which allows general linear constraints (Holmes 2012).  
The vectorized form is:
\deqn{\mathbf{x}(t) = (\mathbf{x}(t-1)^\top \otimes \mathbf{I}_m)(\mathbf{f}_b(t)+\mathbf{D}_b(t)\beta) + (\mathbf{f}_u(t)+\mathbf{D}_u(t)\upsilon) + \mathbf{w}(t), \textrm{ where } \mathbf{W}(t) \sim \textrm{MVN}(0,\textbf{Q}(t))}{x(t) = (t(x(t-1)) (x) I_m)(f_b(t)+D_b(t)beta) + (f_u(t)+D_u(t)upsilon) + w(t), where W(t) ~ MVN(0,Q(t))}
\deqn{\mathbf{y}(t) = (\mathbf{x}(t)^\top \otimes \mathbf{I}_n)(\mathbf{f}_z(t)+\mathbf{D}_z(t)\zeta) + (\mathbf{f}_a(t)+\mathbf{D}_a(t)\alpha) + \mathbf{v}(t), \textrm{ where } \mathbf{V}(t) \sim \textrm{MVN}(0,\textbf{R}(t))}{y(t) = (t(x(t)) (x) I_n)(f_z(t)+D_z(t)zeta) + (f_a(t)+D_a(t)alpha) + v(t), where V(t) ~ MVN(0,R(t))}
\deqn{\mathbf{x}(1) \sim \textrm{MVN}(x0, V0) \textrm{ or } \mathbf{x}(0) \sim \textrm{MVN}(x0, V0)}{x(1) ~ MVN(x0, V0) or x(0) ~ MVN(x0, V0)}
where \eqn{\beta}{beta}, \eqn{\upsilon}{upsilon}, \eqn{\zeta}{zeta}, and \eqn{\alpha}{alpha} are column vectors of estimated values, the \eqn{\mathbf{f}}{f} are column vectors of inputs (fixed values), and the \eqn{\mathbf{D}}{D} are perturbation matrices that align the estimated values into the right rows.  The \eqn{\mathbf{f}}{f} and \eqn{\mathbf{D}}{D} are potentially time-varying. \eqn{\otimes}{(x)} means kronecker product and \eqn{\mathbf{I}_p}{I_p} is a p x p identity matrix.  

Normally the user will specify their model in "marxss" form, perhaps with text short-cuts. The "marxss" form is then converted to "marss" form using the conversion function \code{marxss_to_marss()}. In "marss" form, the D, d, C, and c information is put in A and U respectively. If there are inputs (d and c), then this will make A and U time-varying. This is unfortunate, because this slows down the EM algorithm considerably due to the unfortunate decision (early on) to store time-varying parameters as 3-dimensional. The functions for the "marss" form (in the file \code{MARSS_marss.R}) convert the "marss" form model into vectorized form and prepares the f (fixed) and D (free) matrices that are at the heart of the model specification.  

Note, "marss" is a model form.  A model form is defined by a collection of form functions discussed in \code{\link{marssMODEL}}.  These functions are not exported to the user, but are called by \code{\link{MARSS}()} using the argument \code{form}.  These internal functions convert the users model list into the vectorized form of a MARSS model and do extensive error-checking.   "marxss" is also a model form and these models are also stored in vectorized form (See examples below).

}
\details{
See Holmes (2012) for a discussion of MARSS models in vectorized form. 
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.
}
\seealso{ \code{\link{marssMODEL}}, \code{\link{MARSS.marss}()}, \code{\link{MARSS.marxss}()} } 
\examples{
dat <- t(harborSealWA)
dat <- dat[2:4, ] 
MLEobj <- MARSS(dat)

# free (D) and fixed (f) matrices
names(MLEobj$model$free)
names(MLEobj$model$fixed)
# In marss form, the D, C, d, and c matrices are found in A and U
# If there are inputs, this makes U time-varying
names(MLEobj$marss$free)
names(MLEobj$marss$fixed)

# par is in marss form so does not have values for D, C, d, or c
names(MLEobj$par)
# if you need the par in marxss form, you can use print
tmp <- print(MLEobj, what="par", form="marxss", silent=TRUE)
names(tmp)

}
\references{ 
Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]
}