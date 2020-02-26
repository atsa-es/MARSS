\name{augment.marssMLE}
\alias{augment.marssMLE}

\title{Return the model and state fitted values, residuals, and residual sigma}
\usage{
augment.marssMLE(x, type = c("ytT", "xtT"),
    interval = c("none", "confidence", "prediction"), 
    conf.level = 0.95, form=attr(x[["model"]], "form")[1])
}
\arguments{
  \item{x}{a \code{\link{marssMLE}} object}
  
  \item{type}{What types of fitted values and residuals to return. ytT (observations) and xtT (states) are the values conditioned on all the data.  Read the details below for xtT.  \code{tidy} would be the more common function for returning xtT (smoothed state) estimates.}
  
  \item{interval}{Type of interval: none, confidence or prediction interval.  If the confidence, approximate intervals from the standard errors of the fitted values is given.}
  
  \item{conf.level}{Confidence level.}
  
  \item{form}{If you want the augment function to use a different augment function than augment_form.  This might be useful if you manually specified a DFA model and want to use augment.dfa for rotating.}

}
\description{
\code{augment.marssMLE} returns a data.frame with fitted values, residuals, and upper and lower confidence intervals (if requested) for the fitted observations or states. \code{augment} is concerned with predictions of the states or observations at time t, i.e. the right part of a MARSS model equation with the error terms left off. The error terms are the residuals (the \eqn{v_t}{v(t)} or \eqn{w_t}{w(t)}). \code{\link{tidy.marssMLE}} is concerned with estimates of values (parameters, states or observations) conditioned on all the data.
}
\value{
If interval = "none", the data frame has the following columns:

\item{.fitted}{Fitted values of observations or states. See details.}
\item{.resids}{Model or states residuals. See details.}
\item{.sigma}{The standard error of the model or state residuals. Intervals for the residuals can be constructed from \code{.sigma} using qnorm(alpha/2)*.sigma + .fitted.}
\item{.std.resid}{Standardized residuals. Used for outlier detection. See \code{\link{residuals.marssMLE}}.}

If interval = "confidence", the following are added to the data frame:

\item{.se.fit}{Standard errors of fitted values}
\item{.conf.low}{Lower confidence level at \code{alpha = 1-conf.level}. The interval is approximated using qnorm(alpha/2)*.se.fit + .fitted. }
\item{.conf.up}{Upper confidence level. The interval is approximated using qnorm(1-alpha/2)*.se.fit + .fitted. }

If interval = "prediction", the following are added to the data frame:

\item{.sd.x or .sd.y}{Standard deviation of new x or y iid values.}
\item{.lwr}{Lower range at \code{alpha = 1-conf.level}. The interval is approximated using qnorm(alpha/2)*.sd + .fitted. }
\item{.upr}{Upper range at \code{level}. The interval is approximated using qnorm(1-alpha/2)*.sd + .fitted. }

The standard deviation is for new x or y, not estimated x nor y used in the model. Do not plot the observed data nor the states estimates on the confidence or prediction intervals. For that you need to use \code{.sigma} and construct intervals as noted above. \code{autoplot()} will show the residuals intervals on the observation plot.

}
\details{

See \code{\link{residuals.marssMLE}} for a discussion of the residuals calculations for MARSS models.
The reported CIs are the approximate CIs computed using the standard deviations: \code{qnorm(alpha/2)*se.fit + fitted}.

\strong{observations (type="ytT")}

This returns a model predicted value of the response (\eqn{\mathbf{y}}{y}) and the difference between the model prediction and the observed data is the residual.  If there is no data point, then the residual is NA.  The standard errors help visualize how well the model fits to the data.  See \code{\link{fitted.marssMLE}} for a discussion of the calculation of the fitted values for the observations (the modeled values).  The standardized residuals can be used for outlier detection.  See \code{\link{residuals.marssMLE}} and the chapter on shock detection in the MARSS User Guide.

In the literature on state-space models, it is very common to use the one-step ahead predicted values of the data.  The fitted values returned by \code{type=ytT} are NOT the one-step ahead values and the residuals are NOT the one-step ahead residuals (called Innovations in the state-space literature).  If you want the one-step ahead fitted values, you can use \code{fitted(x, conditioning="t-1")}.  The innovations are also returned by \code{\link{MARSSkf}} in \code{Innov}.

\strong{states (type="xtT")}

If you want the expected value of the states and an estimate of their standard errors (for confidence intervals), then \code{augment} is not what you want to use.  You want to use \code{\link{tidy.marssMLE}} to return the smoothed estimate of the state. \code{augment(MLEobj, type="xtT")} returns a model prediction of \eqn{\mathbf{x}(t)}{x(t)} given \eqn{\mathbf{x}_{t-1}}{x(t-1)}.  The residuals returned are for \eqn{\mathbf{w}_t}{w(t)}, the difference between the two.  These types of residuals are used for outlier detection or shock detection in the state process.  They are also used for model diagnostics.  See \code{\link{residuals.marssMLE}} and read the references cited.

}


\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
MLEobj <- MARSS(dat, model = list(Z = factor(c("WA", "OR", "OR"))))

library(broom)
library(ggplot2)
theme_set(theme_bw())

# Make a plot of the observations and model fits
d <- augment(MLEobj, interval = "confidence")
ggplot(data = d) +
  geom_line(aes(t, .fitted)) +
  geom_point(aes(t, y)) +
  geom_ribbon(aes(x = t, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.1) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count")

# Make a plot of xtT versus prediction of xt from xtT[t-1]
# This is NOT the estimate of the states with CIs. Use tidy() for that.
d <- augment(MLEobj, type = "xtT")
ggplot(data = d) +
  geom_point(aes(t, xtT)) +
  geom_line(aes(x = t, .fitted)) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("xtT (points) and predition (line)")

}

\note{
Within the base code, a form-specific internal \code{augment} function is called to allow the output to vary based on form: augment_dfa, augment_marss, augment_marxss.
}