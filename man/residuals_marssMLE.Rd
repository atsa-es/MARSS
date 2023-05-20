\name{residuals.marssMLE}
\alias{residuals.marssMLE}

\title{Model and state fitted values, residuals, and residual sigma}
\usage{
\method{residuals}{marssMLE}(object, ..., 
           type=c("tt1", "tT", "tt"),
           standardization=c("Cholesky", "marginal", "Block.Cholesky"),
           form=attr(object[["model"]], "form")[1],
           clean=TRUE)
}
\arguments{
  \item{object}{a \code{\link{marssMLE}} object}
  
  \item{type}{\code{tt1} means innovations residuals. The fitted values are computed conditioned on the data up to \eqn{t-1}. See \code{\link[=fitted.marssMLE]{fitted}()} with \code{type="ytt1"} or \code{type="xtt1"}. \code{tT} means smoothation residuals. the fitted values are computed conditioned on all the data. See \code{\link[=fitted.marssMLE]{fitted}()} with \code{type="ytT"} or \code{type="xtT"}. \code{tt} means contemporaneous residuals. The fitted values are computed conditioned on the data up to \eqn{t}. In MARSS functions, estimates at time \eqn{t} conditioned on data 1 to \eqn{T} are denoted \code{tT}, conditioned on the data from \eqn{t=1} to \eqn{t-1} are denoted \code{tt1} and conditioned on data 1 to \eqn{t} are \code{tt}.}

  \item{standardization}{"Cholesky" means it is standardized by the lower triangle of the Cholesky transformation of the full variance-covariance matrix of the model and state residuals. 
"marginal" means that the residual is standardized by its standard deviation, i.e. the square root of the value on the diagonal of the variance-covariance matrix of the model and state residuals. 
"Block.Cholesky" means the model or state residuals are standardized by the lower triangle of the Cholesky transformation of only their variance-covariance matrix (not the joint model and state variance-covariance matrix). }
 
  \item{form}{ For developers. Can be ignored. If you want the function to use a different function than \code{residuals_form}.  This might be useful if you manually specified a DFA model and want to use \code{residuals_dfa} for rotating.}

  \item{clean}{ Can be ignored. For \code{type="tt1"}, state residuals are not used for residuals analysis and for \code{type="tt"}, they don't exist (all NA). They are used only for smoothation residuals, \code{type="tT"}.  For \code{type="tt1"} and \code{type="tt"}, the data frame is cleaned by removing \code{name=="state"} when \code{clean=TRUE}. }
  
  \item{...}{Not used.}

}
\description{
\code{residuals.marssMLE} returns a data frame with fitted values, residuals, residual standard deviation (sigma), and standardized residuals. A residual is the difference between the "value" of the model (\eqn{\mathbf{y}}{y}) or  state (\eqn{\mathbf{x}}{x}) and the fitted value.  At time \eqn{t} (in the returned data frame), the model residuals are for time \eqn{t}. For the the state residuals, the residual is for the transition from \eqn{t} to \eqn{t+1} following the convention in Harvey, Koopman and Penzer (1998). For the the state innovation residuals, this means that \code{state.residual[,t]} is for the transition from \eqn{t} to \eqn{t+1} and is conditioned on data 1 to \eqn{t} while \code{model.residual[,t]} is is conditioned on data 1 to \eqn{t-1}. State innovation residuals are not normally used while state smoothation residuals are used in trend outlier analysis. If warnings are reported, use \code{attr(residuals(fit), "msg")} to retrieve the messages.

Because the state residuals is for the transition from \eqn{t} to \eqn{t+1}, this means that the state residual \code{.resids[t]} is \code{value[t-1]} minus \code{.fitted[t-1]} in the outputted data frame.
}
\value{
A data frame with the following columns:

\item{type}{tT, tt1 or tt}
\item{.rownames}{The names of the observation rows or the state rows.}
\item{name}{model or state}
\item{t}{time step}
\item{value}{The data value if \code{name} equals "model" or the \eqn{x} estimate if \code{name} equals "state" at time \eqn{t}. See details.}
\item{.fitted}{Model predicted values of observations or states at time \eqn{t}. See details.}
\item{.resids}{Model or states residuals. See details.}
\item{.sigma}{The standard error of the model or state residuals. Intervals for the residuals can be constructed from \code{.sigma} using \code{qnorm(alpha/2)*.sigma + .fitted.}}
\item{.std.resids}{Standardized residuals. See \code{\link{MARSSresiduals}} for a discussion of residual standardization.}

}
\details{

See \code{\link{MARSSresiduals}} for a discussion of the residuals calculations for MARSS models.

\strong{value} and \strong{.fitted}

See the discussion below on the meaning of these for \eqn{\mathbf{y}}{y} associated residuals (model residuals) or \eqn{\mathbf{x}}{x} associated residuals (state residuals).

\strong{model residuals}

The model residuals are in the data frame with \code{name=="model"}.

The model residuals are the familiar type of residuals, they are the difference between the data at time \eqn{t} and the predicted value at time \eqn{t}, labeled \code{.fitted} in the data frame. For the model residuals, the "value"" is the data (or NA if data are missing). If \code{type="tT"}, the predicted value is the expected value of \eqn{\mathbf{Y}}{Y} conditioned on all the data, i.e. is computed using the smoothed estimate of \eqn{\mathbf{x}}{x} at time \eqn{t} (\code{xtT}).  If \code{type="tt1"}, the predicted value is the expected value of \eqn{\mathbf{Y}}{Y} conditioned on the data up to time \eqn{t-1}, i.e. is computed using the estimate of \eqn{\mathbf{x}}{x} at time \eqn{t} conditioned on the data up to time \eqn{t-1} (\code{xtt1}). These are known as the one-step-ahead predictions and the residuals are known as the innovations. 

The standard errors help visualize how well the model fits to the data.  See \code{\link[=fitted.marssMLE]{fitted}} for a discussion of the calculation of the model predictions for the observations.  The standardized smoothation residuals can be used for outlier detection.  See the references in \code{\link{MARSSresiduals}} and the chapter on shock detection in the MARSS User Guide.

\strong{state residuals}

The state residuals are in the data frame with \code{name=="state"}.

If you want the expected value of the states and an estimate of their standard errors (for confidence intervals), then \code{residuals()} is not what you want to use.  You want to use \code{\link[=tsSmooth.marssMLE]{tsSmooth}(..., type="xtT")} to return the smoothed estimate of the state or you can find the states in the \code{states} element of the \code{\link{marssMLE}} object returned by a \code{MARSS()} call. For the one-step-ahead state estimates, use \code{\link[=tsSmooth.marssMLE]{tsSmooth}(..., type="xtt1")}.

The state residuals are only for state-space models. At time \eqn{t}, the state residuals are the difference between the state estimate at time \eqn{t+1} and the predicted value of the state at time \eqn{t+1} given the estimate of the state at time \eqn{t}. For smoothation state residuals, this is 
\deqn{ \hat{\mathbf{w}}_{t+1} = \mathbf{x}_{t+1}^T - \mathbf{B}\mathbf{x}_{t}^T - \mathbf{u} - \mathbf{C}\mathbf{c}_{t+1}}{ hatw(t+1) = xtT(t+1) - B xtT(t) - u - C c(t+1)}
For "tt1" state residuals, this is 
\deqn{ \hat{\mathbf{w}}_{t+1} = \mathbf{x}_{t+1}^{t+1} - \mathbf{B}\mathbf{x}_{t}^t - \mathbf{u}  - \mathbf{C}\mathbf{c}_{t+1}}{ hatw(t+1) = xtt(t+1) - B xtt(t) - u - C c(t+1) }.
Note the t indexing is offset. The state residual at time t is the estimate at time t+1 minus the fitted value at t+1.

Smoothation state residuals are used for outlier detection or shock detection in the state process.  See \code{\link{MARSSresiduals}} and read the references cited. Note that the state residual at time \eqn{T} (the last time step) is NA since this would be the transition from \eqn{T} to \eqn{T+1} (past the end of the data).

Note, because the state residuals are for the transition from \eqn{t} to \eqn{t+1}, this means that in the outputted data frame, the state residual \code{.resids[t]} is \code{value[t-1]} minus \code{.fitted[t-1]}.

}


\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
fit <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))

library(ggplot2)
theme_set(theme_bw())

\dontrun{
# Show a series of standard residuals diagnostic plots for state-space models
autoplot(fit, plot.type="residuals")
}

d <- residuals(fit, type="tt1")
\dontrun{
# Make a series of diagnostic plots from a residuals object
autoplot(d)
}

# Manually make a plot of the model residuals (innovations) with intervals
d$.conf.low <- d$.fitted+qnorm(0.05/2)*d$.sigma
d$.conf.up <- d$.fitted-qnorm(0.05/2)*d$.sigma
ggplot(data = d) +
  geom_line(aes(t, .fitted)) +
  geom_point(aes(t, value), na.rm=TRUE) +
  geom_ribbon(aes(x = t, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.1) +
  ggtitle("Model residuals (innovations)") +
  xlab("Time Step") + ylab("Count") +
  facet_grid(~.rownames)

# NOTE state residuals are for t to t+1 while the value and fitted columns 
# are for t. So (value-fitted)[t] matches .resids[t+1] NOT .resids[t]
# This is only for state residuals. For model residuals, the time-indexing matches.
d <- residuals(fit, type="tT")
dsub <- subset(d, name=="state")
# note t in col 1 matches t+1 in col 2
head(cbind(.resids=dsub$.resids, valminusfitted=dsub$value-dsub$.fitted))

# Make a plot of the smoothation residuals
ggplot(data = d) +
  geom_point(aes(t, value-.fitted), na.rm=TRUE) +
  facet_grid(~.rownames+name) +
  ggtitle("Smoothation residuals (state and model)") +
  xlab("Time Step") + ylab("Count")

# Make a plot of xtT versus prediction of xt from xtT[t-1]
# This is NOT the estimate of the smoothed states with CIs. Use tsSmooth() for that.
ggplot(data = subset(d, name=="state")) +
  geom_point(aes(t, value), na.rm=TRUE) +
  geom_line(aes(x = t, .fitted), color="blue") +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("xtT (points) and prediction (line)")

# Make a plot of y versus prediction of yt from xtT[t]
# Why doesn't the OR line go through the points?
# Because there is only one OR state line and it needs to go through
# both sets of OR data.
ggplot(data = subset(d, name=="model")) +
  geom_point(aes(t, value), na.rm=TRUE) +
  geom_line(aes(x = t, .fitted), color="blue") +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("data (points) and prediction (line)")

}

\references{
Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045. 

See also the discussion and references in \code{\link{MARSSresiduals.tT}}, \code{\link{MARSSresiduals.tt1}} and \code{\link{MARSSresiduals.tt}}.
}
