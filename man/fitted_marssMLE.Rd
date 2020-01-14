\name{fitted.marssMLE}
\alias{fitted.marssMLE}

\title{ fitted function for MARSS MLE objects }
\description{
  The MARSS fitting function, \code{\link{MARSS}}, outputs marssMLE objects.  \code{fitted(marssMLE)}, where marssMLE is the output from a \code{\link{MARSS}} call,  will return the modeled value of \eqn{y(t)} or \eqn{x(t)}. For \eqn{y(t)}, this is \eqn{Z x(t) + a}. For \eqn{x(t)}, this is \eqn{B x(t-1) + u}. If you want the estimate of \eqn{x(t)}, then use \code{\link{tidy.marssMLE}}.
}
  \usage{
  \method{fitted}{marssMLE}(object, ..., 
  type=c("observations", "states"),   
  conditioning = c("T", "t", "t-1"),
  interval = c("none", "confidence", "prediction"), 
  conf.level = 0.95,
  output = c("tibble", "matrix"))
  }
  \arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{...}{ Other arguments. Used for backwards compatibility with old arguments. }
  \item{type}{ Fitted values for the observations (y) or the states (x).}
  \item{conditioning}{ The data are t=1 to t=T. If \code{conditioning="T"}, then condition on all the data. If "t", condition on data up to "t". If "t-1", condition on data up to "t-1". The latter is termed "one-step-ahead". Note, for \code{type="states"}, conditioning on data up to "t" is not provided (since this would require the estimate of \eqn{X(t-1)} conditioned on data 1:t).}
  \item{interval}{ If interval="confidence", then the standard error of the fitted value is returned. }
  \item{conf.level}{ Level for the intervals if interval != "none".}
  \item{output}{ tibble or list of matrices}
  }
\value{
  If interval="none" (the default), a T column matrix of the fitted values with one row for each observation (or state) time series is returned. 

If interval = "confidence", the following are returned in a list:

\item{.fitted}{Fitted values of observations or states. See details.}
\item{.se.fit}{Standard errors of fitted values}
\item{.conf.low}{Lower confidence level at \code{alpha = 1-conf.level}. The interval is approximated using qnorm(alpha/2)*.se.fit + .fitted. }
\item{.conf.up}{Upper confidence level. The interval is approximated using qnorm(1-alpha/2)*.se.fit + .fitted. }

The confidence interval is for the fitted value, i.e. the expected value of y (Zx(t)+a) or x (Bx(t-1)+u).

If interval = "prediction", the following are returned in a list:

\item{.fitted}{Fitted values of observations or states. See details.}
\item{.sd.x or .sd.y}{Standard deviation of new x or y iid values.}
\item{.lwr}{Lower range at \code{alpha = 1-conf.level}. The interval is approximated using qnorm(alpha/2)*.sd + .fitted. }
\item{.upr}{Upper range at \code{level}. The interval is approximated using qnorm(1-alpha/2)*.sd + .fitted. }

The prediction interval is for new x or y, not estimated x nor y used in the model. Do not plot the observed data nor the states estimates on these intervals. For that you need the residuals intervals provided by \code{\link{augment.marssMLE}}.

}
\details{
\strong{observation fitted values}

The model predicted (fitted) y at time t is \eqn{Z(t)\hat{x}(t)+a(t)}{Z(t)\hat{x}(t)+a(t)}, where the model is written in marss form. See \code{\link{MARSS.marxss}} for a discussion of the conversion of MARSS models with covariates (c and d) into marss form (which is how models are written in the internal MARSS algorithms).
  
\eqn{\hat{x}(t)} is the expected value of the states at time t. If \code{conditioning="T"}, \eqn{\hat{x}(t)} is the expected value conditioned on all the data, i.e. \code{xtT} returned by \code{\link{MARSSkf}}.  If \code{conditioning="t-1"}, then expected value uses only the data up to time t-1, i.e. \code{xtt1} returned by \code{\link{MARSSkf}}.  These are commonly known as the one step ahead predictions for a state-space model. If \code{conditioning="t"}, then the expected value uses the data up to time t, i.e. \code{xtt} returned by \code{\link{MARSSkf}}. 

If interval="confidence", the se and interval is for the fitted Y. The standard error of the fitted values is Z \%*\% Vt \%*\% t(Z).  If interval="prediction", the sd of new iid Y datasets is returned. The standard deviation of new Y is Z \%*\% Vt \%*\% t(Z) + R.  Vt is either conditioned on 1:T, 1:t, or 1:t-1 depending on \code{conditioning}. Do not plot the data used in the model (to get the estimates of X) on these intervals. For that you want the conditional model residuals (and se's). Use \code{\link{augment.marssMLE}}. See the vignette on residuals in MARSS models. 

\strong{state fitted values}

The model predicted \eqn{\hat{x}(t)} given \eqn{\hat{x}(t-1)} is \eqn{B(t)\hat{x}(t-1)+u(t)}, where the model is written in marss form. This type of state fitted value is used for process outlier detection and shock detection.  See \code{\link{residuals.marssMLE}} and read the references cited.

If you want estimates of the states, rather than the model predictions based on \eqn{x_{t-1}} then you'll want either the states estimate conditioned on all the data (or conditioned on the data up to time t-1 or up to time t).  These are returned by \code{\link{MARSSkf}} in \code{xtT}, \code{xtt1} and \code{xtt} respectively.  Which one you want depends on your objective and application.  You can also use the \code{\link{tidy.marssMLE}} function to return a data.frame with the states (conditioned on data up to T, t, or t-1) with standard errors and intervals.
  
\eqn{\hat{x}(t-1)} used in the prediction is the expected value of the states at time t-1. If \code{conditioning="T"}, \eqn{\hat{x}(t-1)} is the expected value at time t-1 conditioned on all the data, i.e. \code{xtT[,t-1]} returned by \code{\link{MARSSkf}}.  If \code{conditioning="t-1"}, \eqn{\hat{x}(t-1)} is the expected value conditioned on the data up to time t-1, i.e. \code{xtt[,t-1]} returned by \code{\link{MARSSkf}}. The fitted state values conditioned on data up to t is not provided. This would require the expected value of states at time t conditioned on data up to time t+1, and this is not output by the Kalman filter. Only conditioning on data up to t-1 and T are output.

The intervals returned by \code{fitted.marssMLE} for the fitted states, are not typically what one uses or needs--however might be useful for simulation work. If you are doing outlier detection or shock detection, you need the intervals on the state residuals (X(t)-fitted(X(t)) which are returned by \code{\link{residuals.marssMLE}} or \code{\link{augment.marssMLE}}. If you want intervals on the states estimates (conditioned all the data or data up to t-1), use \code{\link{tidy.marssMLE}}.

If interval="confidence", the standard error of the fitted values (meaning the 'expected value of X(t)' given X(t-1)) is returned.  The standard error of the fitted value is B \%*\% Vt[t-1] \%*\% t(B). If interval="prediction", the standard deviation of X(t) given X(t-1) given. The latter is B \%*\% Vt[t-1] \%*\% t(B) + Q (notice it includes Q). Vt is either conditioned on 1:T or 1:t-1 depending on \code{conditioning.} Do not plot (or compare) the estimate of X(t) to the intervals for the fitted values (expected value of X(t) conditioned on X(t-1))! You need the conditional states residuals intervals in that case; use \code{\link{augment.marssMLE}}.
}
\seealso{
\code{\link{MARSSkf}}, \code{\link{augment.marssMLE}}, \code{\link{tidy.marssMLE}}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
MLEobj <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))
fitted(MLEobj)

# Example of fitting a stochastic level model to the Nile River flow data
# red line is smooothed level estimate
# grey line is one-step-ahead prediction using xtt1
nile <- as.vector(datasets::Nile)
mod.list <- list(
  Z = matrix(1), A = matrix(0), R = matrix("r"),
  B = matrix(1), U = matrix(0), Q = matrix("q"),
  x0 = matrix("pi")
)
fit <- MARSS(nile, model = mod.list, silent = TRUE)
plot(nile, type = "p", pch = 16, col = "blue")
lines(fitted(fit)[1, ], col = "red", lwd = 2)
lines(fitted(fit, conditioning = "t-1")[1, ], col = "grey", lwd = 2)
  }