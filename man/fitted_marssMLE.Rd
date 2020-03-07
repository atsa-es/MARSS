\name{fitted.marssMLE}
\alias{fitted.marssMLE}

\title{ fitted function for MARSS MLE objects }
\description{
  \code{\link{MARSS}()} outputs \code{\link{marssMLE}} objects.  \code{fitted(MLEobj)}, where MLEobj is the output from a \code{\link{MARSS}()} call,  will return the modeled value of \eqn{\mathbf{y}_t}{y(t)} or \eqn{\mathbf{x}_t}{x(t)}. For \eqn{\mathbf{y}_t}{y(t)}, this is \eqn{\mathbf{Z}_t\tilde{\mathbf{x}}_t+\mathbf{a}_t}{Z(t) tildex(t)+a(t)}. For \eqn{\mathbf{x}_t}{x(t)}, this is \eqn{\mathbf{B}_t\tilde{\mathbf{x}}_{t-1}+\mathbf{u}_t}{B(t) tildex(t-1)+u(t)}. If you want the estimate of \eqn{\mathbf{x}_t}{x(t)}, then use \code{\link{tidy.marssMLE}}.
}
\usage{
\method{fitted}{marssMLE}(object, ..., 
    type = c("ytT", "xtT", "ytt", "ytt1", "xtt1"),   
    interval = c("none", "confidence", "prediction"), 
    conf.level = 0.95,
    output = c("tibble", "matrix"))
  }
  \arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{...}{ Other arguments. Used for backwards compatibility with old arguments. }
  \item{type}{ Fitted values for the observations (y) or the states (x). If 'tT', then the estimate at time 't' is conditioned on all the data. If 'tt', then the estimate is conditioned on data up to time t. If 'tt1', estimate is conditioned on data up to time t-1. This is also known as one-step-ahead estimate or for y, the innovations.}
  \item{interval}{ If interval="confidence", then the standard error and confidence interval of the fitted value is returned. }
  \item{conf.level}{ Level for the intervals if interval != "none".}
  \item{output}{ tibble or list of matrices}
  }
\value{
  If interval="none" (the default), a T column matrix of the fitted values with one row for each observation (or state) time series is returned. 

If interval = "confidence", the following are returned in a list:

\item{.fitted}{Fitted values of observations (y) or states (x). See details.}
\item{.se.fit}{Standard errors of fitted values}
\item{.conf.low}{Lower confidence level at \code{alpha = 1-conf.level}. The interval is approximated using qnorm(alpha/2)*.se.fit + .fitted. }
\item{.conf.up}{Upper confidence level. The interval is approximated using qnorm(1-alpha/2)*.se.fit + .fitted. }

The confidence interval is for the fitted value, i.e. \eqn{\mathbf{Z}_t\tilde{\mathbf{x}}_t+\mathbf{a}_t}{Z(t) tildex(t)+a(t)} or \eqn{\mathbf{B}_t\tilde{\mathbf{x}}_{t-1}+\mathbf{u}_t}{B(t) tildex(t-1)+u(t)}.

If interval = "prediction", the following are returned in a list:

\item{.fitted}{Fitted values of observations (y) or states (x). See details.}
\item{.sd.x or .sd.y}{Standard deviation of new \eqn{\mathbf{x}_t}{x(t)} or \eqn{\mathbf{y}_t}{y(t)} iid values.}
\item{.lwr}{Lower range at \code{alpha = 1-conf.level}. The interval is approximated using qnorm(alpha/2)*.sd + .fitted. }
\item{.upr}{Upper range at \code{level}. The interval is approximated using qnorm(1-alpha/2)*.sd + .fitted. }

The prediction interval is for new \eqn{\mathbf{x}_t}{x(t)} or \eqn{\mathbf{y}_t}{y(t)}, not estimated \eqn{\mathbf{x}_t}{x(t)} nor \eqn{\mathbf{y}_t}{y(t)} used in the model. Do not plot the observed data nor the states estimates on these intervals. For that you need the residuals intervals provided by \code{\link{augment.marssMLE}}.

}
\details{
\strong{observation fitted values}

The model predicted (fitted) \eqn{\hat{\mathbf{y}}_t}{y(t)} is \eqn{\mathbf{Z}_t\tilde{\mathbf{x}}_t+\mathbf{a}_t}{Z(t) tildex(t)+a(t)}, where the model is written in marss form. See \code{\link{MARSS.marss}} for a discussion of the conversion of MARSS models with covariates (\eqn{\mathbf{c}}{c} and \eqn{\mathbf{d}}{d}) into marss form which is how models are written in the internal MARSS algorithms).
  
\eqn{\tilde{\mathbf{x}}_t}{tildex(t)} is the expected value of the states at time \eqn{t}. If \code{type="ytT"}, \eqn{\tilde{\mathbf{x}}_t}{tildex(t)} is the expected value conditioned on all the data, i.e. \code{xtT} returned by \code{\link{MARSSkf}()}.  If \code{type="ytt1"}, then expected value uses only the data up to time \eqn{t-1}, i.e. \code{xtt1} returned by \code{\link{MARSSkf}()}.  These are commonly known as the one step ahead predictions for a state-space model. If \code{type="ytt"}, then the expected value uses the data up to time \eqn{t}, i.e. \code{xtt} returned by \code{\link{MARSSkf}()}. 

If interval="confidence", the se and interval is for the fitted \eqn{\mathbf{y}}{y}. The standard error of the fitted values is \eqn{\mathbf{Z}_t \tilde{\mathbf{V}}_t \mathbf{Z}_t^\top}{Z(t) tildeV(t) t(Z(t))}.  If interval="prediction", the standard deviation of new iid \eqn{\mathbf{y}}{y} datasets is returned. The standard deviation of new \eqn{\mathbf{y}}{y} is \eqn{\mathbf{Z}_t \tilde{\mathbf{V}}_t \mathbf{Z}_t^\top}{Z(t) tildeVt t(Z(t))}.  \eqn{\tilde{\mathbf{V}}_t}{tildeV(t)} is either conditioned on 1:T, 1:t, or 1:t-1 depending on \code{type}. Prediction intervals are not for the data used in the model but rather new data sets. To evaluate the data used to fit the model, you want the model residuals (and se's). Use \code{\link{augment.marssMLE}} for residuals.

\strong{state fitted values}

The model predicted \eqn{\mathbf{x}(t)}{x(t)} given \eqn{\mathbf{x}_{t-1}}{x(t-1)} is \eqn{\mathbf{B}(t)\tilde{\mathbf{x}}(t-1)+\mathbf{u}(t)}{B(t) tildex(t-1)+u(t)}, where the model is written in "marss" form (\code{\link{MARSS.marss}}). This type of state fitted value is used for process outlier detection and shock detection.  See \code{\link{residuals.marssMLE}} and read the references cited.

If you want estimates of the states, rather than the model predictions based on \eqn{\mathbf{x}_{t-1}}{x(t-1)} then you'll want either the states estimate conditioned on all the data (or conditioned on the data up to time \eqn{t-1} or up to time \eqn{t}).  These are returned by \code{\link{MARSSkf}()} in \code{xtT}, \code{xtt1} and \code{xtt} respectively.  Which one you want depends on your objective and application.  You can also use the \code{\link{tidy.marssMLE}()} function to return a data.frame (tibble) with the estimated states with standard errors and intervals.
  
\eqn{\tilde{\mathbf{x}}_{t-1}}{tildex(t)} used in the prediction is the expected value of the states at time \eqn{t-1}. If \code{type="xtT"}, this is the expected value at time \eqn{t-1} conditioned on all the data, i.e. \code{xtT[,t-1]} returned by \code{\link{MARSSkf}()}.  If \code{type="xtt1"}, it is the expected value conditioned on the data up to time \eqn{t-1}, i.e. \code{xtt[,t-1]} returned by \code{\link{MARSSkf}()}. The fitted state values conditioned on data up to \eqn{t} is not provided. This would require the expected value of states at time \eqn{t} conditioned on data up to time \eqn{t+1}, and this is not output by the Kalman filter. Only conditioning on data up to \eqn{t-1} and T are output.

The intervals returned by \code{fitted.marssMLE} for the fitted states, are not typically what one uses or needs--however might be useful for simulation work. If you are doing outlier detection or shock detection, you need the intervals on the smoothed state residuals which are returned by \code{\link{residuals.marssMLE}} or \code{\link{augment.marssMLE}}. If you want intervals on the states estimates, use \code{\link{tidy.marssMLE}}.

If interval="confidence", the standard error of the fitted values (meaning the 'expected value of \eqn{\mathbf{X}_t}{X(t)}' given \eqn{\mathbf{X}_{t-1}}{X(t-1)}) is returned.  The standard error of the fitted value is \eqn{\mathbf{B}_t \tilde{\mathbf{V}}_{t-1} \mathbf{B}_t^\top}{B(t) tildeV(t-1) t(B(t))}. If interval="prediction", the standard deviation of \eqn{\mathbf{X}_t}{X(t)} given \eqn{\mathbf{X}_{t-1}}{X(t-1)} is output. The latter is \eqn{\mathbf{B}_t \tilde{\mathbf{V}}_{t-1} \mathbf{B}_t^\top + \mathbf{Q}}{B(t) tildeV(t-1) t(B(t)) + Q} (notice it includes Q). \eqn{\tilde{\mathbf{V}}_{t-1}}{tildeV(t-1)} is either conditioned on 1:T or 1:t-1 depending on \code{type}. Do not plot (or compare) the estimate of \eqn{\mathbf{x}_t}{x(t)} to the intervals for the fitted values. You need the conditional states residuals intervals in that case; use \code{\link{augment.marssMLE}}.
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
fit <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))
fitted(fit)

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
lines(fitted(fit, type="ytT")[1, ], col = "red", lwd = 2)
lines(fitted(fit, type="ytt1")[1, ], col = "grey", lwd = 2)
  }