\name{MARSSpredict}
\alias{MARSSpredict}

\title{ Model predictions for MARSS MLE objects }
\description{
This function returns the model predictions for the observations (\eqn{\mathbf{y}}{y}) or the states (\eqn{\mathbf{x}}{x}) conditioned on data from 1 to \eqn{t-1}, \eqn{t} or \eqn{T}. The normal user function for this is the \code{predict} method for  \code{\link{marssMLE}} objects: \code{\link{predict.marssMLE}}.

\code{\link{MARSS}()} outputs \code{\link{marssMLE}} objects.  \code{MARSSpredict(MLEobj)}, where MLEobj is the output from a \code{\link{MARSS}()} call,  will return the modeled value of \eqn{\mathbf{y}_t}{y(t)} or \eqn{\mathbf{x}_t}{x(t)}. For \eqn{\mathbf{y}_t}{y(t)}, this is \eqn{\mathbf{Z}_t\tilde{\mathbf{x}}_t+\mathbf{a}_t}{Z(t) tildex(t)+a(t)}. For \eqn{\mathbf{x}_t}{x(t)}, this is \eqn{\mathbf{B}_t\tilde{\mathbf{x}}_{t-1}+\mathbf{u}_t}{B(t) tildex(t-1)+u(t)}. 

If you want the estimate of \eqn{\mathbf{x}_t}{x(t)} conditioned on all the data or data up to time \eqn{t}, then use \code{\link{fitted.marssMLE}}. Note that the prediction of \eqn{\mathbf{x}_t}{x(t)} conditioned on the data up to time \eqn{t} is not provided since that would require the estimate of \eqn{\mathbf{x}_t}{x(t)} conditioned on data 1 to \eqn{t+1}, which is not output from the Kalman filter or smoother (only conditioning on 1 to \eqn{t-1}, \eqn{t} or \eqn{T} is output).
}
\usage{
MARSSpredict(object, 
    type = c("ytT", "xtT", "ytt", "ytt1", "xtt1"),   
    interval = c("none", "confidence", "prediction"), 
    level = 0.95,
    output = c("data.frame", "matrix"))
  }
  \arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{type}{ If \code{type="tT"}, then the predictions are conditioned on all the data. If  \code{type="tt"}, then the predictions are conditioned on data up to time \eqn{t}. If  \code{type="tt1"}, the predictions are conditioned on data up to time \eqn{t-1}. The latter are also known as one-step-ahead estimate or for \eqn{y}, the innovations.}
  \item{interval}{ If \code{interval="confidence"}, then the standard error and confidence interval of the predicted value is returned. If \code{interval="prediction"}, then the standard deviation and prediction interval of new data or states are returned. }
  \item{level}{ Level for the intervals if interval != "none".}
  \item{output}{ data frame or list of matrices}
  }
\value{
If \code{output="data.frame"} (the default), a data frame with the following columns is returned. If If \code{output="matrix"}, the columns are returned as matrices in a list.

If \code{interval="none"}, the following are returned:
  
\item{.rownames}{ Names of the data or states. }
\item{t}{ Time step. }
\item{y}{ The data if \code{type} is \code{"ytT"}, \code{"ytt"} or \code{"ytt1"}. }
\item{xtT}{ The estimate of \eqn{\mathbf{x}_t}{x_t} conditioned on all the data if \code{type} is \code{"xtT"}. From \code{\link{fitted.marssMLE}()}. }
\item{xtt}{ The estimate of \eqn{\mathbf{x}_t}{x_t} conditioned on the data 1 to \eqn{t} if \code{type} is \code{"xtt1"}. From \code{\link{fitted.marssMLE}()}. }
\item{.pred}{Predicted values of observations (\eqn{\mathbf{y}}{y}) or the states (\eqn{\mathbf{x}}{x}). See details.}


If interval = "confidence", the following are also returned:

\item{.se}{ Standard errors of the predictions. }
\item{.conf.low}{ Lower confidence level at \code{alpha = 1-level}. The interval is approximated using qnorm(alpha/2)*.se + .pred. }
\item{.conf.up}{ Upper confidence level. The interval is approximated using qnorm(1-alpha/2)*.se + .pred. }

The confidence interval is for the predicted value, i.e. \eqn{\mathbf{Z}_t\tilde{\mathbf{x}}_t+\mathbf{a}_t}{Z(t) tildex(t)+a(t)} for \eqn{\mathbf{y}}{y} or \eqn{\mathbf{B}_t\tilde{\mathbf{x}}_{t-1}+\mathbf{u}_t}{B(t) tildex(t-1)+u(t)} for \eqn{\mathbf{x}}{x}.

If interval = "prediction", the following are also returned:

\item{.sd}{ Standard deviation of new \eqn{\mathbf{x}_t}{x(t)} or \eqn{\mathbf{y}_t}{y(t)} iid values. }
\item{.lwr}{ Lower range at \code{alpha = 1-level}. The interval is approximated using qnorm(alpha/2)*.sd + .pred }
\item{.upr}{ Upper range at \code{level}. The interval is approximated using qnorm(1-alpha/2)*.sd + .pred. }

The prediction interval is for new \eqn{\mathbf{x}_t}{x(t)} or \eqn{\mathbf{y}_t}{y(t)}, not estimated \eqn{\mathbf{x}_t}{x(t)} nor \eqn{\mathbf{y}_t}{y(t)} used in the model. Do not plot the observed data nor the states estimates on these intervals. For that you need the residuals intervals provided by \code{\link{residuals.marssMLE}}.

}
\details{
\strong{observation predicted values}

The model predicted \eqn{\hat{\mathbf{y}}_t}{y(t)} is \eqn{\mathbf{Z}_t\tilde{\mathbf{x}}_t+\mathbf{a}_t}{Z(t) tildex(t)+a(t)}, where the model is written in marss form. See \code{\link{MARSS.marss}} for a discussion of the conversion of MARSS models with covariates (\eqn{\mathbf{c}}{c} and \eqn{\mathbf{d}}{d}) into marss form which is how models are written in the internal MARSS algorithms).
  
\eqn{\tilde{\mathbf{x}}_t}{tildex(t)} is the expected value of the states at time \eqn{t}. If \code{type="ytT"}, \eqn{\tilde{\mathbf{x}}_t}{tildex(t)} is the expected value conditioned on all the data, i.e. \code{xtT} returned by \code{\link{MARSSkf}()}.  If \code{type="ytt1"}, then expected value uses only the data up to time \eqn{t-1}, i.e. \code{xtt1} returned by \code{\link{MARSSkf}()}.  These are commonly known as the one step ahead predictions for a state-space model. If \code{type="ytt"}, then the expected value uses the data up to time \eqn{t}, i.e. \code{xtt} returned by \code{\link{MARSSkf}()}. 

If \code{interval="confidence"}, the se and interval is for the predicted \eqn{\mathbf{y}}{y}. The standard error is \eqn{\mathbf{Z}_t \tilde{\mathbf{V}}_t \mathbf{Z}_t^\top}{Z(t) tildeV(t) t(Z(t))}.  If \code{interval="prediction"}, the standard deviation of new iid \eqn{\mathbf{y}}{y} datasets is returned. The standard deviation of new \eqn{\mathbf{y}}{y} is \eqn{\mathbf{Z}_t \tilde{\mathbf{V}}_t \mathbf{Z}_t^\top + \mathbf{R}_t}{Z(t) tildeVt t(Z(t)) + R(t)}.  \eqn{\tilde{\mathbf{V}}_t}{tildeV(t)} is either conditioned on 1 to \eqn{T}, \eqn{t}, or \eqn{t-1} depending on \code{type}. Prediction intervals are not for the data used in the model but rather new data sets. To evaluate the data used to fit the model, you want the model residuals (and residual se's). Use \code{\link{residuals.marssMLE}} for residuals. The intervals for the model residuals are narrower because the fitted \eqn{\mathbf{y}}{y} was estimated from the model data (so is closer to the data than new independent data will be).

\strong{state fitted values}

The model predicted \eqn{\mathbf{x}(t)}{x(t)} given \eqn{\mathbf{x}_{t-1}}{x(t-1)} is \eqn{\mathbf{B}(t)\tilde{\mathbf{x}}(t-1)+\mathbf{u}(t)}{B(t) tildex(t-1)+u(t)}, where the model is written in marss form (\code{\link{MARSS.marss}}). This type of state predicted value is used for process outlier detection and shock detection.  See \code{\link{MARSSresiduals}} and read the references cited.

If you want estimates of the states, rather than the model predictions based on \eqn{\mathbf{x}_{t-1}}{x(t-1)} then you will want either the states estimate conditioned on all the data, on the data up to time \eqn{t-1} or up to time \eqn{t}.  These are returned by \code{\link{fitted.marssMLE}()} in \code{xtT}, \code{xtt1} and \code{xtt} respectively.  Which one you want depends on your objective and application. 
  
\eqn{\tilde{\mathbf{x}}_{t-1}}{tildex(t-1)} used in the prediction is the expected value of the states at time \eqn{t-1}. If \code{type="xtT"}, this is the expected value at time \eqn{t-1} conditioned on all the data, i.e. \code{xtT[,t-1]} returned by \code{\link{MARSSkf}()}.  If \code{type="xtt1"}, it is the expected value conditioned on the data up to time \eqn{t-1}, i.e. \code{xtt[,t-1]} returned by \code{\link{MARSSkf}()}. The predicted state values conditioned on data up to \eqn{t} is not provided. This would require the expected value of states at time \eqn{t} conditioned on data up to time \eqn{t+1}, and this is not output by the Kalman filter. Only conditioning on data up to \eqn{t-1} and \eqn{T} are output.

The intervals returned by \code{MARSSpredict} are for the states predicted based on the state estimate at \eqn{t-1}. These are not typically what one uses or needs--however might be useful for simulation work. If you are doing outlier detection or shock detection, you need the intervals on the smoothed state residuals which are returned by \code{\link{MARSSresiduals}} or \code{\link{residuals.marssMLE}}. If you want intervals on the estimate of the state at time \eqn{t}, use \code{\link{fitted.marssMLE}}.

If \code{interval="confidence"}, the standard error of the predicted values (meaning the standard error of the expected value of \eqn{\mathbf{X}_t}{X(t)} given \eqn{\mathbf{X}_{t-1}}{X(t-1)}) is returned.  The standard error of the predicted value is \eqn{\mathbf{B}_t \tilde{\mathbf{V}}_{t-1} \mathbf{B}_t^\top}{B(t) tildeV(t-1) t(B(t))}. If \code{interval="prediction"}, the standard deviation of \eqn{\mathbf{X}_t}{X(t)} given \eqn{\mathbf{X}_{t-1}}{X(t-1)} is output. The latter is \eqn{\mathbf{B}_t \tilde{\mathbf{V}}_{t-1} \mathbf{B}_t^\top + \mathbf{Q}}{B(t) tildeV(t-1) t(B(t)) + Q} (notice it includes Q). \eqn{\tilde{\mathbf{V}}_{t-1}}{tildeV(t-1)} is either conditioned on 1 to \eqn{T} or \eqn{t-1} depending on \code{type}. Do not plot (or compare) the estimate of \eqn{\mathbf{x}_t}{x(t)} (fitted values) to the intervals for the predicted values. You need the conditional state residuals intervals for that comparison; use \code{\link{residuals.marssMLE}}.
}
\seealso{
\code{\link{MARSSkf}}, \code{\link{MARSSresiduals}}, \code{\link{residuals.marssMLE}}, \code{\link{predict.marssMLE}}, \code{\link{fitted.marssMLE}}
}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
}
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
fit <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))
MARSSpredict(fit)

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
lines(MARSSpredict(fit, type="ytT")$.pred, col = "red", lwd = 2)
lines(MARSSpredict(fit, type="ytt1")$.pred, col = "grey", lwd = 2)
  }