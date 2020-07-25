\name{tsSmooth.marssMLE}
\alias{tsSmooth.marssMLE}
\alias{tsSmooth}

\title{ Smoothed and filtered x and y time series }
\description{


\code{tsSmooth.marssMLE} returns the estimated state and observations conditioned on the data. This function will return either the smoothed values (conditioned on all the data) or the filtered values (conditioned on data 1 to \eqn{t} or \eqn{t-1}). This is output from the Kalman filter and smoother \code{\link{MARSSkf}()} for the \eqn{\mathbf{x}}{x} and from the corresponding function \code{\link{MARSShatyt}()} for the \eqn{\mathbf{y}}{y}. 

These are the expected value of the full right side of the MARSS equations with the error terms (expected value of X(t) and Y(t)). Conditioning on data 1 to \eqn{t-1} (one-step ahead), \eqn{t} (contemporaneous), or \eqn{T} (smoothed) is provided. This is in contrast to \code{\link[dest=fitted.marssMLE]{fitted}} which returns the expected value of the right side without the error term, aka model predictions.

In the state-space literature, the \eqn{\mathbf{y}}{y} "estimates" would normally refer to the expected value of the right-side of the \eqn{\mathbf{y}}{y} equation without the error term (i.e. the expected value of \eqn{\mathbf{Z} \mathbf{X}_t + \mathbf{a}}{Z X(t) + a}). That is provided in \code{\link[dest=fitted.marssMLE]{fitted}}. \code{tsSmooth.marssMLE} provides the expected value with the error terms conditioned on the data from 1 to \eqn{t-1}, \eqn{t}, or \eqn{T}. These estimates are used to estimate missing values in the data. If \eqn{\mathbf{y}}{y} is multivariate, some \eqn{y} are missing at time \eqn{t} and some not, and \eqn{\mathbf{R}}{R} is non-diagonal, then the expected value of \eqn{\mathbf{y}_t}{y(t)} from the right-side of the \eqn{\mathbf{y}}{y} without the error terms would be incorrect because it would not take into account the information in the observed data at time \eqn{t} on the missing data at time \eqn{t} (except as it influences \eqn{\mathbf{x}_t}{x(t)}). 

Note, if there are no missing values, the expected value of \eqn{\mathbf{y}_t}{y(t)} (with error terms) conditioned on the data from 1 to \eqn{t} or \eqn{T} is simply \eqn{\mathbf{y}_t}{y(t)}. The expectation is only useful when there are missing values for which an estimate is needed. The expectation of the \eqn{\mathbf{y}}{y} with the error terms is used in the EM algorithm for the general missing values case and the base function is \code{\link{MARSShatyt}}.
}
\usage{
\method{tsSmooth}{marssMLE}(object, 
    type = c("xtT", "xtt", "xtt1", "ytT", "ytt", "ytt1"),
    interval = c("none", "confidence", "prediction"),
    level = 0.95, ...)
}
\arguments{
  \item{object}{A \code{\link{marssMLE}} object.}
  
  \item{type}{ Type of estimates to return. Smoothed states (xtT), one-step-ahead states (xtt1), contemporaneous states (xtt), the model ytT (Z xtT + A), the model ytt (Z xtt + A)the model ytt1 (Z xtt1 + A), the expected value of y conditioned on data 1 to \eqn{t-1} (ytt1), the expected value of y conditioned on data 1 to \eqn{t} (ytt), or the expected value of y conditioned on data 1 to \eqn{T} (ytT). See details.}
  
\item{interval}{ If \code{interval="confidence"}, then the standard error and confidence intervals are returned. There are no prediction intervals for estimated states and observations except for ytT (which is a unusual case.) If you are looking for prediction intervals, then you want \code{\link[dest=fitted.marssMLE]{fitted}} or \code{\link{predict.marssMLE}}.}

  \item{level}{ Confidence level. alpha=1-level }
 
  \item{...}{Optional arguments.  If form="dfa", \code{rotate=TRUE} can be passed in to rotate the trends (only trends not the Z matrix).}

}

\value{
A data frame with the following columns is returned. Values computed from the model are prefaced with ".".

If \code{interval="none"}, the following are returned:
  
\item{.rownames}{ Names of the data or states. }
\item{t}{ Time step. }
\item{y}{ The data if \code{type} is \code{"ytT"}, \code{"ytt"} or \code{"ytt1"}. }
\item{.estimate}{ The estimated values. See details. }

If \code{interval = "confidence"}, the following are also returned:

\item{.se}{ Standard errors of the estimates. }
\item{.conf.low}{ Lower confidence level at \code{alpha = 1-level}. The interval is approximated using qnorm(alpha/2)*se + estimate }
\item{.conf.up}{ Upper confidence level. The interval is approximated using qnorm(1-alpha/2)*se + estimate }

If \code{interval = "prediction"}, the following are also returned:

\item{.sd}{ Standard deviation of new \eqn{\mathbf{y}_t}{y(t)} values. }
\item{.lwr}{ Lower range at \code{alpha = 1-level}. The interval is approximated using qnorm(alpha/2)*sd + estimate }
\item{.upr}{ Upper range at \code{level}. The interval is approximated using qnorm(1-alpha/2)*sd + estimate }

}

\details{

Below, X and Y refers to the random variable and x and y refer to a specific realization from this random variable.
  
\strong{state estimates (x)}

For \code{type="xtT"}, \code{tsSmooth.marssMLE} returns the confidence intervals of the state at time \eqn{t} conditioned on the data from 1 to \eqn{T} using the estimated model parameters as true values. These are the standard intervals that are shown for the estimated states in state-space models. For example see, Shumway and Stoffer (2000), edition 4, Figure 6.4. As such, this is probably what you are looking for if you want to put intervals on the estimated states (the \eqn{\mathbf{x}}{x}). However, these intervals do not include parameter uncertainty. If you want state residuals (for residuals analysis), use \code{\link{MARSSresiduals}} or \code{\link{residuals.marssMLE}}.

\strong{Quantiles} The state \eqn{\mathbf{X}_t}{X(t)} in a MARSS model has a conditional multivariate normal distribution, that can be computed from the model parameters and data. In Holmes (2012, Equation 11) notation, its expected value conditioned on the observed data and the model parameters \eqn{\Theta}{Theta} is \eqn{\tilde{\mathbf{x}}_t}{tildex(t)}. In \code{\link{MARSSkf}}, this is \code{xtT[,t]}. The variance of \eqn{\mathbf{X}_t}{X(t)} conditioned on the observed data and \eqn{\Theta} is \eqn{\tilde{\mathbf{V}}_t}{tildeV(t)} (\code{VtT[,,t]}). Note that \code{VtT[,,t] != B VtT[,,t-1] t(B) + Q}, which you might think by looking at the MARSS equations. That is because the variance of \eqn{\mathbf{W}_t}{W(t)} conditioned on the data (past, current and FUTURE) is not equal to  \eqn{\mathbf{Q}}{Q} (\eqn{\mathbf{Q}}{Q} is the unconditional variance).

\eqn{\tilde{\mathbf{x}}_t}{tildex(t)} (\code{xtT}) is an estimate of \eqn{\mathbf{x}_t}{x(t)} (the true value), and the standard error of that estimate is given by \eqn{\tilde{\mathbf{V}}_t}{tildeV(t)} (\code{VtT[,,t]}). Let \code{se.xt} denote the sqrt of the diagonal of \code{VtT}. The equation for the \eqn{\alpha/2}{alpha/2} confidence interval is (\code{qnorm(alpha/2)*se.xt + xtT}). \eqn{\mathbf{x}_t}{x(t)} is multivariate and this interval is for one of the \eqn{x}'s in isolation. You could compute the m-dimensional confidence region for the multivariate \eqn{\mathbf{x}_t}{x(t)}, also, but \code{tsSmooth.marssMLE} returns the univariate confidence intervals. 

The variance \code{VtT} gives information on the uncertainty of the true location of \eqn{\mathbf{x}_t}{x(t)} conditioned on the observed data. As more data are collected (or added to the analysis), this variance will shrink since the data, especially data at time \eqn{t}, increases the information about the locations of \eqn{\mathbf{x}_t}{x(t)}. This does not affect the estimation of the model parameters, those are fixed (we are assuming), but rather our information about the states at time \eqn{t}.

If you have a DFA model (form='dfa'), you can pass in \code{rotate=TRUE} to return the rotated trends.  If you want the rotated loadings, you will need to compute those yourself:
\preformatted{
dfa <- MARSS(t(harborSealWA[,-1]), model=list(m=2), form="dfa")
Z.est <- coef(dfa, type="matrix")$Z
H.inv <- varimax(coef(dfa, type="matrix")$Z)$rotmat
Z.rot <- Z.est \%*\% H.inv
}

For \code{type="xtt"} and \code{type=="xtt1"}, the calculations and interpretations of the intervals are the same but the conditioning is for data 1 to \eqn{t} or 1 to \eqn{t-1}. 


\strong{observation estimates (y)}

For \code{type="ytT"}, this returns the expected value and standard error of \eqn{\mathbf{Y}_t}{Y_t} (left-hand side of the \eqn{\mathbf{y}}{y} equation) conditioned on \eqn{\mathbf{Y}_t=y_t}{Y_t=y_t}. If you have no missing data, this just returns your data set. But you have missing data, this what you want in order to estimate the values of missing data in your data set. The expected value of \eqn{\mathbf{Y}_t|\mathbf{Y}=\mathbf{y}(1:T)}{Y|Y=y(1:T)} is in \code{ytT} in \code{\link{MARSShatyt}} output and the variance is \code{OtT-tcrossprod(ytT)} from the \code{\link{MARSShatyt}} output. 

The intervals reported by \code{tsSmooth.marssMLE} for the missing values take into account all the information in the data, specifically the correlation with other data at time \eqn{t} if \eqn{\mathbf{R}}{R} is not diagonal. This is what you want to use for interpolating missing data. You do not want to use \code{predict.marssMLE()} as those predictions are for entirely new data sets and thus will ignore relevant information if \eqn{\mathbf{y}_t}{y_t} is multivariate, not all \eqn{\mathbf{y}_t}{y_t} are missing, and the \eqn{\mathbf{R}}{R} matrix is not diagonal.

The standard error and confidence interval for the expected value of the missing data along with the standard deviation and prediction interval for the missing data are reported. The former uses the variance of \eqn{E[\mathbf{Y}_t]}{E[Y(t)]} conditioned on the data while the latter uses variance of \eqn{\mathbf{Y}_t}{Y(t)} conditioned on the data. \code{\link{MARSShatyt}} returns these variances and expected values. See Holmes (2012) for a discussion of the derivation of expectation and variance of \eqn{\mathbf{Y}_t}{Y(t)} conditioned on the observed data (in the section 'Computing the expectations in the update equations').

For \code{type="ytt"} and \code{type="ytt1"}, only the estimates are provided. \code{\link{MARSShatyt}} does not return the necessary variances matrices for the standard errors for these cases.

}
\references{ 
R. H. Shumway and D. S. Stoffer (2000).  Time series analysis and its applications. Edition 4. Springer-Verlag, New York.

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME]
}
\examples{
dat <- t(harborSeal)
dat <- dat[c(2, 11, 12), ]
fit <- MARSS(dat)

# Make a plot of the estimated states
library(ggplot2)
d <- tsSmooth(fit, type = "xtT", interval="confidence")
ggplot(data = d) +
  geom_line(aes(t, .estimate)) +
  geom_ribbon(aes(x = t, ymin = .conf.low, ymax = .conf.up), linetype = 2, alpha = 0.3) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count")

# Make a plot of the estimates for the missing values
library(ggplot2)
d <- tsSmooth(fit, type = "ytT", interval="confidence")
d2 <- tsSmooth(fit, type = "ytT", interval="prediction")
d$.lwr <- d2$.lwr
d$.upr <- d2$.upr
ggplot(data = d) +
  geom_point(aes(t, .estimate)) +
  geom_line(aes(t, .estimate)) +
  geom_point(aes(t, y), color = "blue", na.rm=TRUE) +
  geom_ribbon(aes(x = t, ymin = .conf.low, ymax = .conf.up), alpha = 0.3) +
  geom_line(aes(t, .lwr), linetype = 2) +
  geom_line(aes(t, .upr), linetype = 2) +
  facet_grid(~.rownames) +
  xlab("Time Step") + ylab("Count") +
  ggtitle("Blue=data, Black=estimate, grey=CI, dash=prediction interval")

# Contrast this with the model prediction of y(t), i.e., put a line through the points
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
  ggtitle("Blue=data, Black=estimate, grey=CI, dash=prediction interval")
}