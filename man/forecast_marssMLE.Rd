\name{forecast.marssMLE}
\alias{forecast.marssMLE}

\title{ forecast function for MARSS MLE objects }
\description{
  \code{\link{MARSS}()} outputs \code{\link{marssMLE}} objects.  \code{forecast(object)}, where object is \code{\link{marssMLE}} object, will return the forecasts of \eqn{\mathbf{y}_t}{y(t)} or the smoothed value of \eqn{\mathbf{x}_t}{x(t)} for \code{h} steps past the end of the model data. \code{forecast(object)} returns a \code{\link{marssPredict}} object which can be passed to \code{\link{plot.marssPredict}} for automatic plotting of the forecast. \code{forecast.marssMLE()} is used by \code{\link{predict.marssMLE}()} to generate forecasts.
}
  \usage{
  \method{forecast}{marssMLE}(object, h=1,
                              conf.level = c(0.80, 0.95),
                              type = c("ytT", "xtT"),
                              newdata = list(y=NULL, c=NULL, d=NULL),
                              interval = c("prediction", "confidence", "none"),
                              fun.kf = c("MARSSkfas", "MARSSkfss"),
                              ...)
   }
  \arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{h}{ Number of steps ahead to forecast. If h !=0, then newdata is for the forecast, i.e. for the h time steps after the end of the model data. See details. }
  \item{conf.level}{ Level for the intervals if interval != "none".}
  \item{type}{ Fitted values for the observations (ytT) or the states (xtT). The predictions are conditioned on all the data (if data are supplied in newdata or implicitly supplied by not passing in newdata).}
  \item{newdata}{ An optional list with new y (data), c or d (covariates) to use for the predictions or forecasts. y, c and d must have the same form as used in the \code{\link{MARSS}()} call but the number of time steps can be different. t is used if there is ambuiquity as to which time steps the y, c and d refer to. See examples and details. }
  \item{interval}{ If interval="confidence", then the standard error and confidence interval of the fitted value is returned. interval="prediction" (default) are the prediction intervals (new data). Note, in the context of a MARSS model, only condidence intervals are available for xtT.}
  \item{fun.kf}{ Only if you want to change the default Kalman filter. Can be ignored. }
  \item{...}{ Other arguments. Not used. }
  }
\value{
  A list with the following components:

\item{method}{The method used for fitting, e.g. MARSS kem. }
\item{model}{The \code{\link{marssMLE}} object passed into predict(). }
\item{newdata}{The \code{newdata} list passed into predict(). c and d must be matrices. t must be a vector (same length as the 2nd dimension of y, c, and d). If h!=0, y must be a matrix. If h==0, y can be a matrix or "none" (no data), or "model" (use the model data in object). If the latter, t can be used to specify a subset of the model data. See examples.}
\item{level}{The \code{conf.level} passed into predict().}
\item{pred}{A data frame the predictions or forecasts along with the intervals. In tibble format. }
\item{type}{The \code{type} ("ytT" or "xtT") passed into predict().}
\item{t}{The time steps used in pred (used for plotting).}
\item{h}{The number of forecast time steps (used for plotting).}


}
\details{

The xtT forecast is returned with \code{\link{tidy.marssMLE}}; it is the smoothed x estimate, i.e. the estimate conditioned on all the data. If h !=0, then 'data' is the original data (model data) plus any new data passed in with y in \code{newdata} for the h forecast time steps. Note, typically forecasts would not have data, since they are forecasts, but \code{predict.marssMLE()} allows you to specify data for the forecast if you need to. If the model includes covariates (c or d matrices passed into the \code{model} list in the \code{\link{MARSS}()} call), then c and d must be passed into \code{newdata}. 

The ytT forecast is returned with \code{\link{fitted.marssMLE}}; it is the fitted ytT or, said alternatively, the expected value of NEW data (y) conditioned on the data used for fitting. The data used for fitting is the same at for xtT (above). The ytT forecast is \eqn{Z xtT + A}. 

If the model has time-varying parameters, the value of the parameters at the last time step is used for the forecast.

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
fr <- predict(models[[i]], type="ytT", h=2)
plot(fr)

# forecast and only show last 10 steps of original data
plot(predict(fit, h=10), include=10)
    
}