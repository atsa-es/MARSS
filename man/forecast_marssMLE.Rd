\name{forecast.marssMLE}
\alias{forecast.marssMLE}

\title{ forecast function for marssMLE objects }
\description{
  \code{\link{MARSS}()} outputs \code{\link{marssMLE}} objects.  \code{forecast(object)}, where object is \code{\link{marssMLE}} object, will return the forecasts of \eqn{\mathbf{y}_t}{y(t)} or \eqn{\mathbf{x}_t}{x(t)} for \code{h} steps past the end of the model data. \code{forecast(object)} returns a \code{\link{marssPredict}} object which can be passed to \code{\link{plot.marssPredict}} for automatic plotting of the forecast. \code{forecast.marssMLE()} is used by \code{\link{predict.marssMLE}()} to generate forecasts.
  
This is a method for the generic \code{forecast} function in the generics package. It is written to mimic the behavior and look of the forecast package.
}
\usage{
\method{forecast}{marssMLE}(object, h = 10, level = c(0.8, 0.95), 
     type = c("ytT","xtT", "ytt", "ytt1", "xtt", "xtt1"), 
     newdata = list(y = NULL, c = NULL, d = NULL), 
     interval = c("prediction", "confidence", "none"), 
     fun.kf = c("MARSSkfas", "MARSSkfss"), ...)
   }
\arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{h}{ Number of steps ahead to forecast. \code{newdata} is for the forecast, i.e. for the \eqn{h} time steps after the end of the model data. If there are covariates in the model, \eqn{\mathbf{c}_t}{c(t)} or \eqn{\mathbf{d}_t}{d(t)}, then \code{newdata} is required. See details. }
  \item{level}{ Level for the intervals if \code{interval != "none"}.}
  \item{type}{ The default for observations would be \code{type="ytT"} and for the states would be \code{type="xtT"}, i.e. using all the data. Other possible forecasts are provided for completeness but would in most cases be identical (see details). }
  \item{newdata}{ An optional list with matrices for new covariates \eqn{\mathbf{c}_t}{c(t)} or \eqn{\mathbf{d}_t}{d(t)} to use for the forecasts. \eqn{\mathbf{c}_t}{c(t)} or \eqn{\mathbf{d}_t}{d(t)} must be in the original model and have the same matrix rows and columns as used in the \code{\link{MARSS}()} call but the number of time steps can be different (and should be equal to \code{h}). }
  \item{interval}{ If \code{interval="confidence"}, then the standard error and confidence interval of the expected value of \eqn{\mathbf{y}_t}{y(t)} (\code{type="ytT"}) or \eqn{\mathbf{x}_t}{x(t)} (\code{type="xtT"}) is returned. \code{interval="prediction"} (default) returns prediction intervals which include the uncertainty in the expected value and due to observation error (the \eqn{\mathbf{R}}{R} in the \eqn{\mathbf{y}}{y} equation). Note, in the context of a MARSS model, only confidence intervals are available for the states (the \eqn{\mathbf{x}}{x}). }
  \item{fun.kf}{ Only if you want to change the default Kalman filter. Can be ignored. }
  \item{...}{ Other arguments. Not used. }
  }
\value{
  A list with the following components:

\item{method}{The method used for fitting, e.g. "kem". }
\item{model}{The \code{\link{marssMLE}} object passed into \code{forecast.marssMLE()}. }
\item{newdata}{The \code{newdata} list if passed into \code{forecast.marssMLE()}.}
\item{level}{The confidence \code{level} passed into \code{forecast.marssMLE()}.}
\item{pred}{A data frame the forecasts along with the intervals. }
\item{type}{The \code{type} ("ytT" or "xtT") passed into \code{forecast.marssMLE()}.}
\item{t}{The time steps used to fit the model (used for plotting).}
\item{h}{The number of forecast time steps (used for plotting).}


}
\details{

The \code{type="ytT"} forecast for \eqn{T+i}  is 
\deqn{\mathbf{Z}\mathbf{x}_{T+i}^T + \mathbf{a} + \mathbf{D}\mathbf{d}_{T+i}}{Z xtT(T+i) + a + D d(T+i)}
where \eqn{\mathbf{Z}}{Z}, \eqn{\mathbf{a}}{a} and \eqn{\mathbf{D}}{D} are estimated from the data from \eqn{t=1} to \eqn{T}. If the model includes \eqn{\mathbf{d}_t}{d(t)} then \code{newdata} with \code{d} must be passed in. Either confidence or prediction intervals can be shown. Prediction intervals would be the norm for forecasts and show the intervals for new data which based on the conditional variance of \eqn{\mathbf{Z}\mathbf{X}_{T+i} + \mathbf{V}_{T+i}}{Z X(T+i) + V(T+i)}. Confidence intervals would show the variance of the mean of the new data (such as if you ran a simulation multiple times and recorded only the mean observation time series). It is based on the conditional variance of \eqn{\mathbf{Z}\mathbf{X}_{T+i}}{Z X(T+i)}. The intervals shown are computed with \code{\link[=fitted.marssMLE]{fitted}()}.

The \code{type="xtT"} forecast for \eqn{T+i}  is 
\deqn{\mathbf{B}\mathbf{x}_{T+i-1}^T + \mathbf{u} + \mathbf{C}\mathbf{c}_{T+i}}{B xtT(T+i-1) + u + C c(T+i)}
where \eqn{\mathbf{B}}{B} and \eqn{\mathbf{u}}{u} and \eqn{\mathbf{C}}{C} are estimated from the data from \eqn{t=1} to \eqn{T} (i.e. the estimates in the marssMLE object). If the model includes \eqn{\mathbf{c}_t}{c(t)} then \code{newdata} with \code{c} must be passed in. The only intervals are confidence intervals which based on the conditional variance of \eqn{\mathbf{B}\mathbf{X}_{T+i-1} + \mathbf{W}_{T+i}}{B X(T+i-1) + W(T+i)}. If you pass in data for your forecast time steps, then the forecast will be computed conditioned on the original data plus the data in the forecast period. The intervals shown are computed with the Kalman smoother (or filter if \code{type="xtt"} or \code{type="xtt1"} specified) via \code{\link[=tsSmooth.marssMLE]{tsSmooth}()}.

If the model has time-varying parameters, the parameter estimates at time \eqn{T} will be used for the whole forecast. If new data \code{c} or \code{d} are passed in, it must have \code{h} time steps. 

Note: \code{y} in \code{newdata}. Data along with covariates can be passed into \code{newdata}. In this case, the data in \code{newdata} (\eqn{T+1} to \eqn{T+h}) are conditioned on for the expected value of \eqn{\mathbf{X}_t}{X(t)}
 but parameters used are only estimated using the data in the marssMLE object (\eqn{t=1} to \eqn{T}). If you include data in \code{newdata}, you need to decide how to condition on that
 new data for the forecast. \code{type="ytT"} would mean that the \eqn{t=T+i} forecast is conditioned on all the data, \eqn{t=1} to \eqn{T+h}, \code{type="ytt"} would mean that the
 \eqn{t=T+i} forecast is conditioned on the data, \eqn{t=1} to \eqn{T+i}, and \code{type="ytt1"} would mean that the \eqn{t=T+i} forecast is conditioned on the data, \eqn{t=1} to \eqn{T+i-1}. Because MARSS models can be used in all sorts of systems, the \eqn{\mathbf{y}}{y} part of the MARSS model might not be "data" in the traditional sense. In some cases, one of the \eqn{\mathbf{y}}{y} (in a multivariate model) might be a known deterministic process or it might be a simulated future \eqn{\mathbf{y}}{y} that you want to include. In this case the
 \eqn{\mathbf{y}}{y} rows that are being forecasted are NAs and the \eqn{\mathbf{y}}{y} rows that are known are passed in with \code{newdata}.
}
\seealso{
\code{\link{plot.marssPredict}}, \code{\link{predict.marssMLE}}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\examples{
# More examples are in ?predict.marssMLE

dat <- t(harborSealWA)
dat <- dat[2:4,] #remove the year row
fit <- MARSS(dat, model=list(R="diagonal and equal"))

# 2 steps ahead forecast
fr <- forecast(fit, type="ytT", h=2)
plot(fr)

# forecast and only show last 10 steps of original data
fr <- forecast(fit, h=10)
plot(fr, include=10)

}