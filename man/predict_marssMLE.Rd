\name{predict.marssMLE}
\alias{predict.marssMLE}

\title{ predict and forecast MARSS MLE objects }
\description{
This function will return the modeled value of \eqn{\mathbf{y}_t}{y(t)} or \eqn{\mathbf{x}_t}{x(t)} conditioned on the data (either the data used to fit the model or data in \code{newdata}). For \eqn{\mathbf{y}_t}{y(t)}, this is \eqn{\mathbf{Z}_t \mathbf{x}_t^T+\mathbf{a}_t+\mathbf{D}_t\mathbf{d}_t}{Z(t) xtT(t) + a(t) + D(t) d(t)}. For \eqn{\mathbf{x}_t}{x(t)}, this is \eqn{\mathbf{B}_t \mathbf{x}_{t-1}^T+\mathbf{u}_t+\mathbf{C}_t\mathbf{c}_{t}}{B(t) xtT(t-1) + u(t) + C(t) c(t)}.  \eqn{\mathbf{x}_t^T}{xtT(t)} is the smoothed state estimate at time \eqn{t} conditioned on all the data (either data used to fit the model or the optional data passed into \code{newdata}).

If you want the estimate of \eqn{\mathbf{x}_t}{x(t)} conditioned on all the data (i.e. output from the Kalman filter or smoother), then use \code{\link[=tsSmooth.marssMLE]{tsSmooth}()}. Note that the prediction of \eqn{\mathbf{x}_t}{x(t)} conditioned on the data up to time \eqn{t} is not provided since that would require the estimate of \eqn{\mathbf{x}_t}{x(t)} conditioned on data 1 to \eqn{t+1}, which is not output from the Kalman filter or smoother.

If \code{h} is passed in, \code{predict(object)} will return a forecast \eqn{h} steps past the end of the model data. \code{predict(object)} returns a \code{\link{marssPredict}} object which can be passed to \code{plot()} or \code{ggplot2::autoplot()}for automatic plotting of predictions and forecasts with intervals.
}
\usage{
\method{predict}{marssMLE}(object, n.ahead = 0, 
     level = c(0.80, 0.95), 
     type = c("ytt1", "ytT", "xtT", "ytt", "xtt1"), 
     newdata = list(t=NULL, y=NULL, c=NULL, d=NULL), 
     interval = c("none", "confidence", "prediction"), 
     fun.kf = c("MARSSkfas", "MARSSkfss"), 
     x0 = "reestimate", ...)
   }
\arguments{
\item{object}{ A \code{\link{marssMLE}} object.  }
\item{n.ahead}{ Number of steps ahead to forecast. If \code{n.ahead !=0}, then \code{newdata} is for the forecast, i.e. for the n.ahead time steps after the end of the model data. See details. }
\item{level}{ Level for the intervals if \code{interval != "none".}}
\item{type}{ \code{ytT}, \code{ytt} or \code{ytt1}: predictions for the observations based on the states estimate at time \eqn{t} conditioned on all the data, data up to \eqn{t} or data up to \eqn{t-1}. \code{xtT} or \code{xtt1}: predictions for the states at time \eqn{t} based on the states estimate at time \eqn{t-1} conditioned on all the data or data up to \eqn{t-1}. The data are the data used to fit the model unless \code{y} is passed in in \code{newdata}.}
\item{newdata}{ An optional list with new \code{y} (data), \code{c} or \code{d} (covariates) to use for the predictions or forecasts. \code{y}, \code{c} or \code{d} must have the same structure (matrix dimensions) as used in the \code{\link{MARSS}()} call but the number of time steps can be different. \code{t} is used if there is ambuiquity as to which time steps the \code{newdata} refer to. See examples and details. }
\item{interval}{ If \code{interval="confidence"}, then the standard error and confidence intervals of the predictions are returned. If \code{interval="prediction"}, prediction intervals are returned. See \code{\link[=fitted.marssMLE]{fitted}} for a discussion of the intervals. }
\item{fun.kf}{ Only if you want to change the default Kalman filter. Can be ignored. }
\item{x0}{ If "reestimate" (the default), then the initial value for the states is re-estimated. If \code{"use.model"}, then the initial values in the fitted model (\code{object}) are used. If you change the data, then this initial condition may not be appropriate. You can also pass in a new \code{x0} to use. It must be a matrix that is the same dimensions as \code{x0} in the model. \code{x0} is ignored if \code{h!=0} since in that case a forecast is being done. See example. }
\item{...}{ Other arguments. Not used. }
}

\value{
  A list with the following components:

\item{method}{The method used for fitting, e.g. MARSS kem. }
\item{model}{The \code{\link{marssMLE}} object passed into predict(). }
\item{newdata}{The \code{newdata} list if passed in.}
\item{level}{The confidence or prediction intervals \code{level}.}
\item{pred}{A data frame the predictions or forecasts along with the intervals. }
\item{type}{The \code{type} passed in.}
\item{t}{The time steps in the pred data frame.}
\item{n.ahead and h}{The number of forecast time steps.}
\item{x0}{The x0 used for the predictions.}
\item{tinitx}{The tinitx used.}

The pred data frame has the following columns:

\item{.rownames}{ Names of the data or states. }
\item{t}{ Time step. }
\item{y}{ The data if \code{type} is \code{"ytT"}, \code{"ytt"} or \code{"ytt1"}. }
\item{xtT}{ The estimate of \eqn{\mathbf{x}_t}{x_t} conditioned on all the data if \code{type="xtT"}. From \code{\link[=tsSmooth.marssMLE]{tsSmooth}()}. }
\item{xtt}{ The estimate of \eqn{\mathbf{x}_t}{x_t} conditioned on the data 1 to \eqn{t} if \code{type="xtt1"}. From \code{\link[=tsSmooth.marssMLE]{tsSmooth}()}. }
\item{estimate}{Model predicted values of observations (\eqn{\mathbf{y}}{y}) or the states (\eqn{\mathbf{x}}{x}). See details.}

If intervals are returned, the following are added to the data frame:

\item{se}{ Standard errors of the predictions. }
\item{Lo ...}{ Lower confidence level at \code{alpha = 1-level}. The interval is approximated using qnorm(alpha/2)*se + prediction. }
\item{Hi ...}{ Upper confidence level. The interval is approximated using qnorm(1-alpha/2)*se + prediction. }

}

\details{

\strong{Forecasts \code{n.ahead != 0}}

The \code{type="xtT"} forecast is the states forecast conditioned on all the data. If \code{n.ahead !=0}, then 'data' that is being conditioned on is the original data (model data) plus any data in \code{newdata$y} for the h forecast time steps. Note, typically forecasts would not have data, since they are forecasts, but \code{predict.marssMLE()} allows you to specify data for the forecast time steps if you need to. If the model includes covariates (\eqn{\mathbf{c}}{c} and/or \eqn{\mathbf{d}}{d} matrices passed into the \code{model} list in the \code{\link{MARSS}()} call), then \code{c} and/or \code{d} must be passed into \code{newdata}. 

The \code{type="ytT"} forecast is the expected value of NEW data (\eqn{\mathbf{y}}{y}) conditioned on the data used for fitting. The data used for fitting is the same as for \code{type="xtT"} (above). The \eqn{\mathbf{y}}{y} forecast is \code{Z xtT[,T+i] + A + D d[,T+i]}. 

If the model has time-varying parameters, the value of the parameters at the last time step are used for the forecast.

\strong{Model predictions \code{n.ahead == 0}}

If \code{newdata} is not passed in, then the model data (\eqn{\mathbf{y}}{y}) and \eqn{\mathbf{c}}{c} and \eqn{\mathbf{d}}{d} (if part of model) are used for the predictions. \code{\link[=fitted.marssMLE]{fitted}(object, type="ytT")} is the internal function for model predictions in that case.

If \code{newdata} is passed in, then the predictions are computed using \code{newdata} but with the MARSS model estimated from the original data, essentially the Kalman filter/smoother is run using the estimated MARSS model but with data (and \eqn{\mathbf{c}}{c} and \eqn{\mathbf{d}}{d} if in the model) in \code{newdata}.  \code{y}, \code{c} and \code{d} in the \code{newdata} list must all have the same number of columns (time-steps) and the length of \code{t} in \code{newdata} must be the same as the number of columns and must be sequential.

For \code{type="ytT"}, the predictions are conceptually the same as predictions returned by \code{predict.lm} for a linear regression. The confidence interval is the interval for the expected value of NEW data. The prediction interval is the interval for NEW data. Prediction intervals will always be wider (or equal if R=0) to confidence intervals. The difference is that the uncertainty in \code{predict.lm} comes from parameter uncertainty and the data error while in \code{predict.marssMLE}, the uncertainty is from \eqn{\mathbf{x}}{x} uncertainty and data error. Parameter uncertainty does not enter the interval calculations; parameters are treated as known at their point estimates. This is not specific to the MARSS package. This is how prediction and confidence intervals are presented for MARSS models in the literature, i.e. no parameter uncertainty. 

\itemize{
\item{\code{t} in \code{newdata}:}{ If the model has time-varying parameters, \code{t} in \code{newdata} removes any ambiguity as to which parameter values (time steps) will be used for prediction. In this case, \code{t} specifies which time values of the parameters you want to use. If you leave off \code{t}, then it is assumed that \code{t} starts at the first time step in the data used to fit the original model. If the model is time-constant, \code{t} is used to set the time step values (used for plotting, etc.). }
\item{The model has \eqn{\mathbf{c}}{c} and/or \eqn{\mathbf{d}}{d}:}{ \code{c} and/or \code{d} must be included in \code{newdata}. If \code{y} (new data) is not in \code{newdata}, it is assumed to be absent (all NA). That is, the default behavior if \code{y} is absent but \code{c} and/or \code{d} is present is \code{y="none"}. If you want to use the original data used to fit the model, then pass in \code{y="model"} in \code{newdata}. Pass in \code{t} in \code{newdata} if it is ambiguous which time steps of the model data to use.}
\item{The model has time-varying parameters:}{ You have to pass in \code{t} in \code{newdata} to specify what parameter values to use. If any \eqn{t > T} (\eqn{T} equals the last time step in the model data), then it is assumed that you want to use the parameter values at the last time step of the original time series for values beyond the last time step. See examples. }
\item{\code{y}, \code{c} and \code{d} in \code{newdata} have more time steps than the original data: }{If the model has time-varying parameters, you will need to pass in \code{t}. If the model is time-constant, then \code{t} is assumed to start at the first time step in the original data but you can pass in \code{t} to change that. It will not change the prediction, but will change the t column in the output.}
}

\strong{x0 estimation} If you are passing in \code{y} in \code{newdata}, then it is likely that you will need to re-estimate the \eqn{\mathbf{x}}{x} initial condition. The default behavior of \code{predict.marssMLE}. Use \code{x0 = "use.model"} to use the initial values in the estimated model (\code{object}).


}
\seealso{
\code{\link{plot.marssPredict}}, \code{\link{fitted.marssMLE}}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
}
\examples{
dat <- t(harborSealWA)
dat <- dat[2:4,] #remove the year row
fit <- MARSS(dat, model=list(R="diagonal and equal"))

# 2 steps ahead forecast
fr <- predict(fit, type="ytT", n.ahead=2)
plot(fr)

# use model data with the estimated initial values (at t=0) for 
# initial values at t=9
# This would be a somewhat strange thing to do and the value at t=10 will look wrong.
fr <- predict(fit,  newdata=list(t=10:20, y=dat[,10:20]), x0 = "use.model")
plot(fr)

# pass in new data and give it new t; initial conditions will be estimated
fr <- predict(fit,  newdata=list(t=23:33, y=matrix(10,3,11)))
plot(fr, ylim=c(8,12))


# Covariate example
fulldat <- lakeWAplanktonTrans
years <- fulldat[,"Year"]>=1965 & fulldat[,"Year"]<1975
dat <- t(fulldat[years,c("Greens", "Bluegreens")])
dat <- zscore(dat)
covariates <- rbind(
  Temp = fulldat[years, "Temp"],
  TP = fulldat[years, "TP"])
covariates <- zscore(covariates)
A <- U <- "zero"
B <- Z <- "identity"
R <- diag(0.16,2)
Q <- "equalvarcov"
C <- "unconstrained"
model.list <- list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=covariates)
fit <- MARSS(dat, model=model.list)

# Use a new c (covariate) but no data.
fr <- predict(fit,  newdata=list(c=matrix(5,2,10)), x0="use.model")
plot(fr)

# Use first 10 time steps of model data
plot(predict(fit, newdata=list(y=dat[,1:10], c=matrix(5,2,10))))

# Use all model data but new covariates
# Why does it look so awful? Because this is a one-step ahead
# prediction and there is no info on what the c will be at t
plot(predict(fit,  newdata=list(y=dat, c=matrix(5,2,120))))

# Use all model data but new covariates with ytT type
# this looks better because is uses all the c data to estimate (so knows what c is at t and beyond)
plot(predict(fit,  newdata=list(y=dat, c=matrix(5,2,120)), type="ytT"))

# Use no data; cannot estimate initial conditions without data
# so x0 must be "use.model"
fr <- predict(fit,  newdata=list(c=matrix(5,2,22)), x0="use.model")
plot(fr)

# forecast with covariates
# n.ahead and the number column in your covariates in newdata must match
plot(predict(fit,  newdata=list(c=matrix(5,2,10)), n.ahead=10))

# forecast with covariates and only show last 10 steps of original data
plot(predict(fit,  newdata=list(c=matrix(5,2,10)), n.ahead=10), include=10)
    
}