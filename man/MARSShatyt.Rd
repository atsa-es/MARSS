\name{MARSShatyt}
\alias{MARSShatyt}
\title{ Compute Expected Value of Y, YY, and YX}
\description{
Computes the expected value of random variables involving Y for the EM algorithm.  Users can also use \code{print( MLEobj, what="Ey")} to access this output.  See \code{\link{print.marssMLE}}.
}
\usage{
MARSShatyt( MLEobj )
}
\arguments{
  \item{ MLEobj }{ A \code{\link{marssMLE}} object with the \code{par} element of estimated parameters, \code{model} element with the model description and data. }
  }
\details{
For state space models, \code{MARSShatyt()} computes the expectations involving Y.  If Y is completely observed, this entails simply replacing Y with the observed y.  When Y is only partially observed, the expectation involves the conditional expectation of a multivariate normal. 
 
}
\value{ 
A list with the following components (n is the number of state processes). Following the notation in Holmes (2012), y(1) is the observed data (for t=1:TT) while y(2) is the unobserved data.  y(1,1:t) is the observed data from time 1 to t.
  \item{ytT}{ Estimates E[Y(t) | Y(1,1:TT)=y(1,1:TT)] (n x T matrix). }
  \item{ytt1}{ Estimates E[Y(t) | Y(1,1:t-1)=y(1,1:t-1)] (n x T matrix). }
  \item{OtT}{ Estimates E[Y(t) t(Y(t) | Y(1,1:TT)=y(1,1:TT)] (n x n x T array). }
  \item{Ott1}{ Estimates E[Y(t) t(Y(t) | Y(1,1:t-1)=y(1,1:t-1)] (n x n x T array). }
  \item{yxtT}{ Estimates E[Y(t) t(X(t) | Y(1,1:TT)=y(1,1:TT)] (n x m x T array). }
  \item{yxt1T}{ Estimates E[Y(t) t(X(t-1) | Y(1,1:TT)=y(1,1:TT)] (n x m x T array). }
  \item{yxttpT}{ Estimates E[Y(t) t(X(t+1) | Y(1,1:TT)=y(1,1:TT)] (n x m x T array). }
  \item{errors}{ Any error messages due to ill-conditioned matrices. }
  \item{ok}{ (TRUE/FALSE) Whether errors were generated. }
}
\references{
Holmes, E. E. (2012) Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models. Technical report. arXiv:1302.3919 [stat.ME] Type \code{RShowDoc("EMDerivation",package="MARSS")} to open a copy.
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ 
\code{\link{MARSS}} \code{\link{marssMODEL}} \code{\link{MARSSkem}}
}


