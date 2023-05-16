\name{MARSSinits}
\alias{MARSSinits}
\keyword{helper}
\title{ Initial Values for MLE }
\description{
  Sets up generic starting values for parameters for maximum-likelihood estimation algorithms that use an iterative maximization routine needing starting values.   Examples of such algorithms are the EM algorithm in \code{\link{MARSSkem}()} and Newton methods in \code{\link{MARSSoptim}()}.   This is a utility function in the \code{\link{MARSS-package}}.  It is not exported to the user. Users looking for information on specifying initial conditions should look at the help file for \code{\link{MARSS}()} and the User Guide section on initial conditions.
  
The function assumes that the user passed in the inits list using the parameter names in whatever form was specified in the \code{\link{MARSS}()} call.  The default is form="marxss".  The \code{\link{MARSSinits}()} function calls MARSSinits_foo, where foo is the form specified in the \code{\link{MARSS}()} call.  MARSSinits_foo translates the inits list in form foo into form marss.
}
\usage{
MARSSinits(MLEobj, inits=list(B=1, U=0, Q=0.05, Z=1, A=0, 
    R=0.05, x0=-99, V0=5, G=0, H=0, L=0))
}
\arguments{
  \item{MLEobj}{ An object of class \code{\link{marssMLE}}. }
  \item{inits}{ A list of column vectors (matrices with one column) of the estimated values in each parameter matrix.}
}
\details{
  Creates an \code{inits} parameter list for use by iterative maximization algorithms.  
  
  Default values for \code{inits} is supplied in \code{MARSSsettings.R}. The user can alter these and supply any of the following (m is the dim of X and n is the dim of Y in the MARSS model): 
\itemize{
    \item{elem=\code{A,U}} A numeric vector or matrix which will be constructed into \code{inits$elem} by the command \code{array(inits$elem),dim=c(n or m,1))}. If elem is fixed in the model, any \code{inits$elem} values will be overridden and replaced with the fixed value. Default is \code{array(0,dim=c(n or m,1))}. 
    \item{elem=\code{Q,R,B}} A numeric vector or matrix.  If length equals the length \code{MODELobj$fixed$elem} then \code{inits$elem} will be constructed by \code{array(inits$elem),dim=dim(MODELobj$fixed$elem))}. If length is 1 or equals dim of \code{Q} or dim of \code{R} then \code{inits$elem} will be constructed into a diagonal matrix by the command \code{diag(inits$elem)}. If elem is fixed in the model, any \code{inits$elem} values will be overridden and replaced with the fixed value. Default is \code{diag(0.05, dim of Q or R)} for \code{Q} and \code{R}.  Default is \code{diag(1,m)} for \code{B}.
    \item{\code{x0}} If \code{inits$x0=-99}, then starting values for \code{x0} are estimated by a linear regression 
    through the count data assuming \code{A} is all zero.  This will be a poor start if \code{inits$A} is not 0. If \code{inits$x0} is a numeric vector or matrix, \code{inits$x0} will be constructed by the command \code{array(inits$x0),dim=c(m,1))}. If \code{x0} is fixed in the model, any \code{inits$x0} values will be overridden and replaced with the fixed value. Default is \code{inits$x0=-99}. 
    \item{\code{Z}} If \code{Z} is fixed in the model, \code{inits$Z} set to the fixed value. If \code{Z} is not fixed, then the user must supply \code{inits$Z}.  There is no default.
    \item{elem=\code{V0}} \code{V0} is never estimated, so this is never used.
}
}
\value{
  A list with initial values for the estimated values for each parameter matrix in a MARSS model in marss form. So this will be a list with elements \code{B}, \code{U}, \code{Q}, \code{Z}, \code{A}, \code{R}, \code{x0}, \code{V0}, \code{G}, \code{H}, \code{L}.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\seealso{ \code{\link{marssMODEL}}, \code{\link{MARSSkem}()}, \code{\link{MARSSoptim}()} }

\note{
Within the base code, a form-specific internal \code{MARSSinits} function is called to allow the output to vary based on form: \code{MARSSinits_dfa}, \code{MARSSinits_marss}, \code{MARSSinits_marxss}.
}
