\name{MARSShatyt}
\alias{MARSShatyt}
\title{ Compute Expected Value of Y, YY, and YX}
\description{
Computes the expected value of random variables involving Y.  Users can also use \code{print( MLEobj, what="Ey")} to access this output.  See \code{\link{print.marssMLE}}.
}
\usage{
MARSShatyt(MLEobj, only.kem = TRUE)
}
\arguments{
  \item{ MLEobj }{ A \code{\link{marssMLE}} object with the \code{par} element of estimated parameters, \code{model} element with the model description and data. }
  \item{ only.kem }{ Return only ytT, OtT, yxtT, and yxttpT (values conditioned on the data from 1:T). If \code{only.kem} If TRUE, only return (and compute) values needed for the EM algorithm.  If only.kem=FALSE, then also return values conditioned on data from 1 to t-1 (Ott1 and ytt1) and 1 to t (Ott and ytt), yxtt1T (E[Y(t),X(t-1)|1:T]), var.ytT (var[Y(t)|1:T]), and var.EytT (var_X[E_{Y|x}[Y(t)|1:T,x(t)]]).}
  }
\details{
For state space models, \code{MARSShatyt()} computes the expectations involving Y.  If Y is completely observed, this entails simply replacing Y with the observed y.  When Y is only partially observed, the expectation involves the conditional expectation of a multivariate normal. 
 
}
\value{ 
A list with the following components (n is the number of state processes). Following the notation in Holmes (2012), y(1) is the observed data (for t=1:T) while y(2) is the unobserved data.  y(1,1:t-1) is the observed data from time 1 to t-1.
  \item{ytT}{ E[Y(t) | Y(1,1:T)=y(1,1:T)] (n x T matrix). }
  \item{ytt1}{ E[Y(t) | Y(1,1:t-1)=y(1,1:t-1)] (n x T matrix). }
  \item{ytt}{ E[Y(t) | Y(1,1:t)=y(1,1:t)] (n x T matrix). }
  \item{OtT}{ E[Y(t) t(Y(t)) | Y(1,1:T)=y(1,1:T)] (n x n x T array). }
  \item{var.ytT}{ var[Y(t) | Y(1,1:T)=y(1,1:T)] (n x n x T array). }
  \item{var.EytT}{ var_X[E_{Y|x}[Y(t) | Y(1,1:T)=y(1,1:T), X(t)=x(t)]] (n x n x T array). }
  \item{Ott1}{  E[Y(t) t(Y(t)) | Y(1,1:t-1)=y(1,1:t-1)] (n x n x T array). }
  \item{Ott}{ E[Y(t) t(Y(t)) | Y(1,1:t)=y(1,1:t)] (n x n x T array). }
  \item{yxtT}{ E[Y(t) t(X(t)) | Y(1,1:T)=y(1,1:T)] (n x m x T array). }
  \item{yxtt1T}{ E[Y(t) t(X(t-1)) | Y(1,1:T)=y(1,1:T)] (n x m x T array). }
  \item{yxttpT}{ E[Y(t) t(X(t+1)) | Y(1,1:T)=y(1,1:T)] (n x m x T array). }
  \item{errors}{ Any error messages due to ill-conditioned matrices. }
  \item{ok}{ (TRUE/FALSE) Whether errors were generated. }
}
\references{
Holmes, E. E. (2012) Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models. Technical report. arXiv:1302.3919 [stat.ME] Type \code{RShowDoc("EMDerivation",package="MARSS")} to open a copy. See the section on 'Computing the expectations in the update equations' and the subsections on expectations involving Y. 
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ 
\code{\link{MARSS}}, \code{\link{marssMODEL}}, \code{\link{MARSSkem}}
}
\examples{
dat <- t(harborSeal)
dat <- dat[2:3, ]
MLEobj <- MARSS(dat)
EyList <- MARSShatyt(MLEobj)
}


