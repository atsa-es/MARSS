\name{forecast.marssMLE}
\alias{forecast.marssMLE}

\title{ forecast function for MARSS MLE objects }
\description{
  \code{\link{MARSS}()} outputs \code{\link{marssMLE}} objects.  \code{forecast(object)}, where object is \code{\link{marssMLE}} object, will return the forecasts of \eqn{\mathbf{y}_t}{y(t)} or \eqn{\mathbf{x}_t}{x(t)} for \code{h} steps past the end of the model data. \code{forecast(object)} returns a \code{\link{marssPredict}} object which can be passed to \code{\link{plot.marssPredict}} for automatic plotting of the forecast. \code{forecast.marssMLE()} is used by \code{\link{predict.marssMLE}()} to generate forecasts.
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
  \item{h}{ Number of steps ahead to forecast. \code{newdata} is for the forecast, i.e. for the \eqn{h} time steps after the end of the model data. If there are covariates, \eqn{c} or \eqn{d}, then \code{newdata} is required. See details. }
  \item{level}{ Level for the intervals if \code{interval != "none"}.}
  \item{type}{ Observations (ytT) or the states (xtT). Other possible forecasts are provided for completeness. }
  \item{newdata}{ An optional list with matrices for new covariates \eqn{c} or \eqn{d} to use for the forecasts. \eqn{c} and \eqn{d} must be in the original model and have the same matrix rows and columns as used in the \code{\link{MARSS}()} call but the number of time steps can be different. }
  \item{interval}{ If \code{interval="confidence"}, then the standard error and confidence interval of the expected value of \eqn{\mathbf{y}_t}{y(t)} (ytT) or \eqn{\mathbf{x}_t}{x(t)} (xtT) is returned. \code{interval="prediction"} (default) returns prediction intervals which include the uncertainty in the expected value and due to observation error (the \eqn{R} in the \eqn{y} equation). Note, in the context of a MARSS model, only confidence intervals are available for x. }
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

The ytT forecast for \eqn{T+i}  is 
\deqn{Z\tilde{x}_{T+i} + a + Dd_{T+i}}
where \eqn{Z}, \eqn{a} and \eqn{D} are estimated form the data from 1 to \eqn{T}. If the model does include \eqn{d} then \code{newdata} with \eqn{d} must be passed in. Either confidence or prediction intervals can be shown. Prediction intervals would be the norm for forecasts and show the intervals for new data which based on the conditional variance of \eqn{ZX_{T+i} + V_{t+i}}. Confidence intervals would show the variance of the mean of the new data (such as if you ran a simulation multiple times and recoded only the mean observation time series). It is based on the conditional variance of \eqn{ZX_{T+i}}. The intervals shown are computed with \code{\link[=fitted.marssMLE]{fitted}()}.

The xtT forecast for \eqn{T+i}  is 
\deqn{B\tilde{x}_{T+i-1} + u + Cc_{T+i}}
where \eqn{B} and \eqn{u} and \eqn{C} are estimated form the data from 1 to \eqn{T} (i.e. the estimates in the marssMLE object). \eqn{Cc_{T+i}} is only present if the model includes \eqn{c}. If the model does include \eqn{c} then \code{newdata} with \eqn{c} must be passed in. The only intervals are confidence intervals which based on the conditional variance of \eqn{BX_{T+i-1} + W_{t+i}}. If you pass in data for your forecast time steps, then the forecast will be computed conditioned on the original data plus the data in the forecast period. The intervals shown are computed with the Kalman smoother (or filter if \code{type="xtt"} or \code{type="xtt1"} specified) via \code{\link[=tsSmooth.marssMLE]{tsSmooth}()}.

If the model has time-varying parameters, the parameter estimates at time \eqn{T} will be used for the whole forecast. If new data \eqn{c} or \eqn{d} is passed in, it must have \eqn{h} time steps. 

Note: \eqn{y} in \code{newdata}. Data along with covariates can be passed into \code{newdata}. In this case, the data in \code{newdata} (\eqn{T+1} to \eqn{T+h}) are conditioned on for the expected value of \eqn{\mathbf{x}_t}{x(t)} 
 but parameters used are only estimated using the data in the marssMLE object (t=1 to T). The expected value of \eqn{y} conditioned on observed \eqn{y} is simply \eqn{y}.
 
}
\seealso{
\code{\link{plot.marssPredict}}, \code{\link{predict.marssMLE}()}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{
# More examples are in ?predict.marssMLE

dat <- t(harborSealWA)
dat <- dat[2:4,] #remove the year row
fit <- MARSS(dat, model=list(R="diagonal and equal"))

# 2 steps ahead forecast
fr <- forecast.marssMLE(fit, type="ytT", h=2)
plot(fr)

# forecast and only show last 10 steps of original data
fr <- forecast.marssMLE(fit, h=10)
plot(fr, include=10)

# if you load the forecast package, you can use
library(forecast)
fr <- forecast(fit)

}