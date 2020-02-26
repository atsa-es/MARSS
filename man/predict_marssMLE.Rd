\name{predict.marssMLE}
\alias{predict.marssMLE}

\title{ predict function for MARSS MLE objects }
\description{
  \code{\link{MARSS}()} outputs \code{\link{marssMLE}} objects.  \code{predict(object)}, where object is \code{\link{marssMLE}} object,  will return the fitted value of \eqn{\mathbf{y}_t}{y(t)} or the smoothed value of \eqn{\mathbf{x}_t}{x(t)}. If \code{h} is passed in, \code{predict(object)} will return a forecast h steps past the end of the model data. \code{predict(object)} returns a \code{\link{marssPredict}} object which can be passed to \code{plot()} for automatic plotting of predictions and forecasts with intervals.
}
  \usage{
  \method{predict}{marssMLE}(object, h=0,
                             conf.level = c(0.80, 0.95),
                             type = c("ytT", "xtT"),
                             newdata = list(t=NULL, y=NULL, c=NULL, d=NULL),
                             interval = c("prediction", "confidence", "none"),
                             fun.kf = c("MARSSkfas", "MARSSkfss"),
                             x0 = list(x0=NULL, tinitx=NULL),
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
  \item{x0}{ If new data (y) is passed in in newdata, you may want to re-estimate the initial condition. If so, pass in the x0 specification in the x0 list. See example. }
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

\strong{Forecasts h != 0}

The xtT forecast is returned with \code{\link{tidy.marssMLE}}; it is the smoothed x estimate, i.e. the estimate conditioned on all the data. If h !=0, then 'data' is the original data (model data) plus any new data passed in with y in \code{newdata} for the h forecast time steps. Note, typically forecasts would not have data, since they are forecasts, but \code{predict.marssMLE()} allows you to specify data for the forecast if you need to. If the model includes covariates (c or d matrices passed into the \code{model} list in the \code{\link{MARSS}()} call), then c and d must be passed into \code{newdata}. 

The ytT forecast is returned with \code{\link{fitted.marssMLE}}; it is the fitted ytT or, said alternatively, the expected value of NEW data (y) conditioned on the data used for fitting. The data used for fitting is the same at for xtT (above). The ytT forecast is \eqn{Z xtT + A}. 

If the model has time-varying parameters, the value of the parameters at the last time step is used for the forecast.

\strong{Predictions h == 0}

If newdata is not passed in, then the model data (y) and c and d (if part of model) are used for the predictions. The xtT is returned with \code{\link{tidy.marssMLE}} and ytT, with prediction or confidence intervals, is returned with \code{\link{fitted.marssMLE}}.

If newdata is passed in, then xtT and ytT are computed using newdata but with the MARSS model estimated from the original data, essentially the Kalman filter/smoother is run using the MARSS model in object but with data (and c and d if in the model) in newdata.  y, c and d must all have the same number of columns (time-steps) and the length of t must be the same as the number of columns and must be sequential.

\itemize{
\item{t in newdata}{t in newdata removes any ambiguity as to which time steps your y, c and d correspond to. In most cases when you pass in y, c or d for a prediction, you will want to specify t (and in some cases need to specify t). If you pass in y and all the parameters are time-constant, t is only used for plotting (the x axis). If you pass in y="model", then t can be used to specify a subset of the model data. If left off, it is assumed you want all the model data. If the model has time-varying parameters, t specifies which values of the parameters you want. If you leave off t, then it is assumed that t is 1 to the number of columns in y, c and d.}
\item{The model has c or d}{c or d must be passed in. If y (new data) is not passed in, it is assumed to be absent. That is, the default behavior if y is absent but c or d is present is y="none". If you want to use the model data, then pass in y="model" in newdata. Pass in t in newdata if it is ambiguous which time steps of the model data to use.}
\item{The model has time-varying parameters}{You have to pass in t in newdata to specify what parameter values to use. Any t > TT (time steps in the model data), is assumed to mean that you want to use the parameter values at the at the last time step of the original time series. See examples. }
\item{Only t is passed in}{You can specify a subset of the model data (y) (and c and d if in model), by passing in t alone with the desired time steps denoted. For example t=1:10.}
\item{y, c and d have more time steps than the original data}{If the model has time-varying parameters, you will need to pass in t. If the model is time-constant, then t is assumed to be 1 to the number of columns in y, c and d but you can pass in t to change that.}
}

\strong{x0 estimation} If you are using different data, then it is likely that you will need to re-estimate x0. You can specify how to estimate x0 as you would in the \code{\link{MARSS}()} \code{model} argument. See examples.


}
\seealso{
\code{\link{plot.marssPredict}}, \code{\link{fitted.marssMLE}}, \code{\link{tidy.marssMLE}}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{
dat <- t(harborSealWA)
dat <- dat[2:4,] #remove the year row
fit <- MARSS(dat, model=list(R="diagonal and equal"))

# 2 steps ahead forecast
fr <- predict(fit, type="ytT", h=2)
plot(fr)

# use model data with initial values (at t=0) for initial values at t=9
fr <- predict(fit,  newdata=list(t=10:20), use.initial.values = TRUE))
plot(fr)

# pass in new data and give it new t; initial conditions will be estimated
fr <- predict(fit,  newdata=list(t=23:33, y=matrix(10,3,11)))
plot(fr, ylim=c(8,12))


# Covariate example
fulldat = lakeWAplanktonTrans
years = fulldat[,"Year"]>=1965 & fulldat[,"Year"]<1975
dat = t(fulldat[years,c("Greens", "Bluegreens")])
dat = zscore(dat)
covariates = rbind(
  Temp = fulldat[years,"Temp"],
  TP = fulldat[years,"TP"])
covariates = zscore(covariates)
A = U = "zero"
B = Z = "identity"
R = diag(0.16,2)
Q = "equalvarcov"
C = "unconstrained"
model.list = list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=covariates)
fit = MARSS(dat, model=model.list)

# Use first 10 time steps of model data
plot(predict(fit, newdata=list(t=1:10, y="model", c=matrix(5,2,10))))

# Use all model data
plot(predict(fit,  newdata=list(t=1:22, y="model", c=matrix(5,2,22))))

# Use no data; cannot estimate initial conditions without data
fr <- predict(fit,  newdata=list(y="none", c=matrix(5,2,22)), use.initial.values=TRUE)
plot(fr)

# forecast with covariates
plot(predict(fit,  newdata=list(c=matrix(5,2,10)), h=10))

# forecast with covariates and only show last 10 steps of original data
plot(predict(fit,  newdata=list(c=matrix(5,2,10)), h=10), include=10)
    
}