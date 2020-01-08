\name{MARSS.out}
\alias{MARSS.out}
\title{ Output and Intervals for MARSS MLE objects }
\description{
  There are a variety of types of outputs and intervals for MARSS models. This help file summarizes them and tells you which function to use to get the output you want.
}
\usage{
MARSS.out() 
}
\details{
When looking at estimates and intervals for x and y, keep in mind that the MARSS parameters are always fixed. The uncertainty arise because x and y are stochastic not because of uncertainty in the parameter estimates. The latter type of uncertainty arises because we could have collected different data and those data would give us different parameter estimates. If you looked at the intervals on model fits (a regression line say) from a classic linear regression those intervals reflect parameter uncertainty. The intervals on x (states) and y (observations) in a MARSS model are fundamentally different. X is stochastic and unobserved. Y is stochastic and observed (at least partially).
\itemize{
\item  If you want the estimate of x conditioned on all the data 1 to T, then use \code{\link{tidy.marssMLE}} with \code{type="states", conditioning="T"}. This gives \eqn{E[X|all data]}.  This is the states estimates from the Kalman smoother. You will can also output these \code{\link{MARSSkf}}. The confidence interval is based on  \eqn{\var[X|all data]}, the conditional variance of X. That variance is \eqn{VtT[,,t]}.

\item  If you want the estimate of x conditioned on the data up to t-1, then use \code{\link{tidy.marssMLE}} with \code{type="states", conditioning="t-1"}. This is the one step ahead prediction, \eqn{B x(t-1) + u}, and its variance. \code{\link{tidy.marssMLE}} returns a prediction interval based on the variance  \eqn{Vtt1[,,t]=B Vtt[,,t-1] B^T + Q}. If you wanted the confidence interval on the expected value, that would use the variance \eqn{B Vtt[,,t-1] B^T}. \code{\link{augment.marssMLE}} will return that.

\item  If you want the estimate of y conditioned on all the data, then use \code{\link{tidy.marssMLE}} with \code{type="observations", conditioning="T"}. This gives \eqn{E[Y|all data]} and if you have no missing data, \eqn{E[Y|data]} = \eqn{data} so you will just be returning your data. The confidence intervals will also be 0. For missing y values, you will see prediction intervals since \code{\link{tidy.marssMLE}} returns intervals based on \eqn{var[Y|all data]}.
  
\item  If you want the estimate of y conditioned on the data up to t-1, then use \code{\link{tidy.marssMLE}} with \code{type="observations", conditioning="t-1"}. This is the one step ahead prediction, \eqn{Z x(t|t-1) + A}, and its variance. \code{\link{tidy.marssMLE}} returns a prediction interval  based on the variance \eqn{Z Vtt1[,,t] Z^T + R}. If you wanted the confidence interval on the expected value of the prediction, that would use \eqn{Z Vtt1[,,t] Z^T}. \code{\link{augment.marssMLE}} will return that.
}

}
\seealso{
\code{\link{MARSSkf}}, \code{\link{augment.marssMLE}}, \code{\link{tidy.marssMLE}}, \code{\link{fitted.marssMLE}}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
