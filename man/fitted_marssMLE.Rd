\name{fitted.marssMLE}
\alias{fitted.marssMLE}
\title{Return fitted values for X(t) and Y(t) in a MARSS model}
\description{

\code{fitted()} returns the different types of fitted values for \eqn{\mathbf{x}_t}{x(t)} and \eqn{\mathbf{y}_t}{y(t)} in a MARSS model. The fitted values are the expected value of the right side of the MARSS equations without the error terms, thus are the model predictions of \eqn{\mathbf{y}_t}{y(t)} or \eqn{\mathbf{x}_t}{x(t)}. \code{fitted.marssMLE} is a companion function to \code{\link[=tsSmooth.marssMLE]{tsSmooth}()} which returns the expected value of the right side of the MARSS equations with the error terms (the Kalman filter and smoother output).
\deqn{\mathbf{x}_{t} = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{C} \mathbf{c}_t + \mathbf{G} \mathbf{w}_t, \textrm{ where } \mathbf{w}_t \sim \textrm{MVN}(0,\mathbf{Q})}{x(t) = B x(t-1) + u + C c(t) + G w(t)}
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{D} \mathbf{d}_t + \mathbf{H} \mathbf{v}_t, \textrm{ where } \mathbf{v}_t \sim \textrm{MVN}(0,\mathbf{R})}{y(t) = Z x(t) + a + D d(t) + H v(t)}
The data go from \eqn{t=1} to \eqn{t=T}. For brevity, the parameter matrices are shown without a time subscript, however all parameters can be time-varying.

Note that the prediction of \eqn{\mathbf{x}_t}{x(t)} conditioned on the data up to time \eqn{t} is not provided since that would require the expected value of \eqn{\mathbf{X}_{t}}{X(t)} conditioned on data from \eqn{t = 1} to \eqn{t+1}, which is not output from the Kalman filter or smoother.
}
\usage{
\method{fitted}{marssMLE}(object, ...,
    type = c("ytt1", "ytT", "xtT", "ytt", "xtt1"),   
    interval = c("none", "confidence", "prediction"), 
    level = 0.95,
    output = c("data.frame", "matrix"),
    fun.kf = c("MARSSkfas", "MARSSkfss"))
  }
  \arguments{
  \item{object}{ A \code{\link{marssMLE}} object.  }
  \item{type}{ If \code{type="tT"}, then the predictions are conditioned on all the data. If  \code{type="tt"}, then the predictions are conditioned on data up to time \eqn{t}. If  \code{type="tt1"}, the predictions are conditioned on data up to time \eqn{t-1}. The latter are also known as one-step-ahead estimates. For \eqn{\mathbf{y}}{y}, these are also known as the innovations.}
  \item{interval}{ If \code{interval="confidence"}, then the standard error and confidence interval of the predicted value is returned. If \code{interval="prediction"}, then the standard deviation and prediction interval of new data or states are returned. }
  \item{level}{ Level for the intervals if \code{interval} is not equal to \code{"none"}.}
  \item{output}{ data frame or list of matrices}
  \item{fun.kf}{ By default, \code{tsSmooth()} will use the Kalman filter/smoother function in \code{object$fun.kf} (either \code{\link{MARSSkfas}()} or \code{\link{MARSSkfss}()}). You can pass in \code{fun.kf} to force a particular Kalman filter/smoother function to be used. }
  \item{...}{ Not used. }
  }
\value{
If \code{output="data.frame"} (the default), a data frame with the following columns is returned. If If \code{output="matrix"}, the columns are returned as matrices in a list. Information computed from the model has a leading "." in the column name.

If \code{interval="none"}, the following are returned (colnames with \code{.} in front are computed values):
  
\item{.rownames}{ Names of the data or states. }
\item{t}{ Time step. }
\item{y}{ The data if \code{type} is \code{"ytT"}, \code{"ytt"} or \code{"ytt1"}. }
\item{.x}{ The expected value of \eqn{\mathbf{X}_t}{X_t} conditioned on all the data if \code{type="xtT"} or data up to time \eqn{t} if \code{type="xtt1"}. From \code{\link[=tsSmooth.marssMLE]{tsSmooth}()}. This is the expected value of the right-side of the \eqn{\mathbf{x}_t}{x(t)} equation with the errors terms while \code{.fitted} is the expected value of the right side without the error term \eqn{\mathbf{w}_t}{w(t)}.}
\item{.fitted}{Predicted values of observations (\eqn{\mathbf{y}}{y}) or the states (\eqn{\mathbf{x}}{x}). See details.}


If \code{interval = "confidence"}, the following are also returned:

\item{.se}{ Standard errors of the predictions. }
\item{.conf.low}{ Lower confidence level at \code{alpha = 1-level}. The interval is approximated using qnorm(alpha/2)*.se + .fitted }
\item{.conf.up}{ Upper confidence level. The interval is approximated using qnorm(1-alpha/2)*.se + .fitted }

The confidence interval is for the predicted value, i.e. \eqn{\mathbf{Z}\mathbf{x}_t^\tau+\mathbf{a}}{Z xttau(t)+a} for \eqn{\mathbf{y}}{y} or \eqn{\mathbf{B}\mathbf{x}_{t-1}^\tau+\mathbf{u}}{B xttau(t-1)+u} for \eqn{\mathbf{x}}{x} where \eqn{\mathbf{x}_t^\tau}{xttau(t)} is the expected value of \eqn{\mathbf{X}_t}{X(t)} conditioned on the data from 1 to \eqn{\tau}{tau}. (\eqn{\tau}{tau} will be \eqn{t-1}, \eqn{t} or \eqn{T}).

If \code{interval = "prediction"}, the following are also returned:

\item{.sd}{ Standard deviation of new \eqn{\mathbf{x}_t}{x(t)} or \eqn{\mathbf{y}_t}{y(t)} iid values. }
\item{.lwr}{ Lower range at \code{alpha = 1-level}. The interval is approximated using qnorm(alpha/2)*.sd + .fitted }
\item{.upr}{ Upper range at \code{level}. The interval is approximated using qnorm(1-alpha/2)*.sd + .fitted }

The prediction interval is for new \eqn{\mathbf{x}_t}{x(t)} or \eqn{\mathbf{y}_t}{y(t)}. If you want to evaluate the observed data or the states estimates for outliers then these are not the intervals that you want. For that you need the residuals intervals provided by \code{\link[=residuals.marssMLE]{residuals}()}.

}
\details{

In the state-space literature, the two most commonly used fitted values are \code{"ytt1"} and 
\code{"ytT"}. The former is the expected value of \eqn{\mathbf{Y}_t}{Y(t)} conditioned on the data 1 to time \eqn{t-1}. These are known as the innovations and they, plus their variance, are used in the calculation of the likelihood of a MARSS model via the innovations form of the likelihood. The latter, \code{"ytT"} are the model estimates of the \eqn{\mathbf{y}}{y} values using all the data; this is the expected value of \eqn{\mathbf{Z}\mathbf{X}_t+\mathbf{a}+\mathbf{D}\mathbf{d}_t}{Z X(t) + a + D d(t)} conditioned on the data 1 to \eqn{T}. The \code{"xtT"} along with \code{"ytT"} are used for computing smoothation residuals used in outlier and shock detection. See \code{\link{MARSSresiduals}}. For completeness, \code{fitted.marssMLE} will also return the other possible model predictions with different conditioning on the data (1 to \eqn{t-1}, \eqn{t}, and \eqn{T}), however only \code{type="ytt1"} (innovations) and \code{"ytT"} and \code{"xtT"} (smoothations) are regularly used. Keep in mind that the fitted \code{"xtT"} is not the smoothed estimate of \eqn{\mathbf{x}} (unlike \code{"ytT"}). If you want the smoothed estimate of \eqn{\mathbf{x}} (i.e. the expected value of \eqn{\mathbf{X}_t}{X(t)} conditioned on all the data), you want the Kalman smoother. See \code{\link[=tsSmooth.marssMLE]{tsSmooth}}.

Fitted versus estimated values: The fitted states at time \eqn{t} are predictions from the estimated state at time \eqn{t-1} conditioned on the data: expected value of \eqn{\mathbf{B}\mathbf{X}_{t-1}+\mathbf{u}+\mathbf{C}\mathbf{c}_t}{B X(t-1) + u + C c(t)} conditioned on the data. They are distinguished from the estimated states at time \eqn{t} which would includes the conditional expected values of the error terms: \eqn{E[\mathbf{X}_{t}] = \mathbf{B}\mathbf{X}_{t-1}+\mathbf{u}+\mathbf{C}\mathbf{c}_t + \mathbf{w}_t}{E[X(t)] = B X(t-1) + u + C c(t) + w(t)}. The latter are returned by the Kalman filter and smoother. Analogously, the fitted observations at time \eqn{t} are model predictions from the estimated state at time \eqn{t} conditioned on the data: the expected value of the right side of the \eqn{\mathbf{y}}{y} equation without the error term. Like with the states, one can also compute the expected value of the observations at time \eqn{t} conditioned on the data: the expected value of the right side of the \eqn{\mathbf{y}}{y} equation with the error term.  The latter would just be equal to the data if there are no missing data, but when there are missing data, this is what is used to estimate their values. The expected value of states and observations are provided via \code{\link[=tsSmooth.marssMLE]{tsSmooth}}.

\strong{observation fitted values}

The model predicted \eqn{\hat{\mathbf{y}}_t}{y(t)} is \eqn{\mathbf{Z}\mathbf{x}_t^\tau+\mathbf{a} + \mathbf{D}\mathbf{d}_t}{Z xttau(t) + a + D d(t)}, where \eqn{\mathbf{x}_t^\tau}{xttau(t)} is the expected value of the state at time \eqn{t} conditioned on the data from 1 to \eqn{\tau}{tau} (\eqn{\tau}{tau} will be \eqn{t-1}, \eqn{t} or \eqn{T}). Note, if you are using MARSS for estimating the values for missing data, then you want to use \code{\link[=tsSmooth.marssMLE]{tsSmooth}()} with \code{type="ytT"} not \code{fitted(..., type="ytT")}.
  
\eqn{\mathbf{x}_t^\tau}{xttau(t)} is the expected value of the states at time \eqn{t} conditioned on the data from time 1 to \eqn{\tau}{tau}. If \code{type="ytT"}, the expected value is conditioned on all the data, i.e. \code{xtT} returned by \code{\link{MARSSkf}()}.  If \code{type="ytt1"}, then the expected value uses only the data up to time \eqn{t-1}, i.e. \code{xtt1} returned by \code{\link{MARSSkf}()}.  These are commonly known as the one step ahead predictions for a state-space model. If \code{type="ytt"}, then the expected value uses the data up to time \eqn{t}, i.e. \code{xtt} returned by \code{\link{MARSSkf}()}. 

If \code{interval="confidence"}, the standard error and interval is for the predicted \eqn{\mathbf{y}}{y}. The standard error is \eqn{\mathbf{Z} \mathbf{V}_t^\tau \mathbf{Z}^\top}{Z Vttau(t) t(Z)}.  If \code{interval="prediction"}, the standard deviation of new iid \eqn{\mathbf{y}}{y} data sets is returned. The standard deviation of new \eqn{\mathbf{y}}{y} is \eqn{\mathbf{Z} \mathbf{V}_t^\tau \mathbf{Z}^\top + \mathbf{R}_t}{Z Vttau(t) t(Z) + R(t)}.  \eqn{\mathbf{V}_t^\tau}{Vttau(t)} is conditioned on the data from \eqn{t=1} to \eqn{n}. \eqn{\tau}{tau} will be either \eqn{t}, \eqn{t-1} or \eqn{T} depending on the value of \code{type}. 

Intervals returned by \code{MARSSpredict()} are not for the data used in the model but rather new data sets. To evaluate the data used to fit the model for residuals analysis or analysis or model inadequacy, you want the model residuals (and residual se's). Use \code{\link{residuals.marssMLE}} for model residuals and their intervals. The intervals for model residuals are narrower because the predictions for \eqn{\mathbf{y}}{y} were estimated from the model data (so is closer to the data used to estimate the predictions than new independent data will be).

\strong{state fitted values}

The model predicted \eqn{\mathbf{x}_t}{x(t)} given \eqn{\mathbf{x}_{t-1}}{x(t-1)} is \eqn{\mathbf{B}\mathbf{x}_{t-1}^\tau+\mathbf{u}+\mathbf{C}\mathbf{c}_t}{B xttau(t-1) + u + C c(t)}. If you want estimates of the states, rather than the model predictions based on \eqn{\mathbf{x}_{t-1}}{x(t-1)}, go to \code{\link[=tsSmooth.marssMLE]{tsSmooth}()}. Which function you want depends on your objective and application. 
  
\eqn{\mathbf{x}_{t-1}^\tau}{xttau(t-1)} used in the prediction is the expected value of the states at time \eqn{t-1} conditioned on the data from $t=1$ to $t=\tau$. If \code{type="xtT"}, this is the expected value at time \eqn{t-1} conditioned on all the data, i.e. \code{xtT[,t-1]} returned by \code{\link{MARSSkf}()}.  If \code{type="xtt1"}, it is the expected value conditioned on the data up to time \eqn{t-1}, i.e. \code{xtt[,t-1]} returned by \code{\link{MARSSkf}()}. The predicted state values conditioned on data up to \eqn{t} are not provided. This would require the expected value of states at time \eqn{t} conditioned on data up to time \eqn{t+1}, and this is not output by the Kalman filter. Only conditioning on data up to \eqn{t-1} and \eqn{T} are output.

If \code{interval="confidence"}, the standard error of the predicted values (meaning the standard error of the expected value of \eqn{\mathbf{X}_t}{X(t)} given \eqn{\mathbf{X}_{t-1}}{X(t-1)}) is returned.  The standard error of the predicted value is \eqn{\mathbf{B} \mathbf{V}_{t-1}^\tau\mathbf{B}^\top}{B Vttau(t-1) t(B)}. If \code{interval="prediction"}, the standard deviation of \eqn{\mathbf{X}_t}{X(t)} given \eqn{\mathbf{X}_{t-1}}{X(t-1)} is output. The latter is \eqn{\mathbf{B} \mathbf{V}_{t-1}^\tau \mathbf{B}^\top + \mathbf{Q}}{B Vttau(t-1) t(B) + Q} . \eqn{\mathbf{V}_{t-1}^\tau}{Vttau(t-1)} is either conditioned on data 1 to \eqn{\tau=T} or \eqn{\tau=t-1} depending on \code{type}. 

The intervals returned by \code{MARSSpredict()} are for the state predictions based on the state estimate at \eqn{t-1}. These are not typically what one uses or needs--however might be useful for simulation work. If you want confidence intervals for the state estimates, those are provided in \code{\link[=fitted.marssMLE]{fitted}()}. If you want to do residuals analysis (for outliers or model inadequacy), you want the residuals intervals provided in \code{\link[=residuals.marssMLE]{residuals}()}.

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
lines(fitted(fit, type="ytT")$.fitted, col = "red", lwd = 2)
lines(fitted(fit, type="ytt1")$.fitted, col = "grey", lwd = 2)

# Make a plot of the model estimate of y(t), i.e., put a line through the points
# Intervals are for new data not the blue dots 
# (which were used to fit the model so are not new)
library(ggplot2)
d <- fitted(fit, type = "ytT", interval="confidence", level=0.95)
d2 <- fitted(fit, type = "ytT", interval="prediction", level=0.95)
d$.lwr <- d2$.lwr
d$.upr <- d2$.upr
ggplot(data = d) +
  geom_line(aes(t, .fitted), size=1) +
  geom_point(aes(t, y), color = "blue", na.rm=TRUE) +
  geom_ribbon(aes(x = t, ymin = .conf.low, ymax = .conf.up), alpha = 0.3) +
  geom_line(aes(t, .lwr), linetype = 2) +
  geom_line(aes(t, .upr), linetype = 2) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("Blue=data, Black=estimate, grey=CI, dash=prediction interval") +
  geom_text(x=15, y=7, label="The intervals are for \n new data not the blue dots")

}