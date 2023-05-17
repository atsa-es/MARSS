\name{MARSSresiduals.tt1}
\alias{MARSSresiduals.tt1}

\title{ MARSS One-Step-Ahead Residuals }
\description{
  Calculates the standardized (or auxiliary) one-step-ahead residuals, aka the innovations residuals and their variance.  Not exported. Access this function with \code{MARSSresiduals(object, type="tt1")}. To get the residuals as a data frame in long-form, use \code{\link[=residuals.marssMLE]{residuals}(object, type="tt1")}.
}
\usage{
MARSSresiduals.tt1(object, method = c("SS"), normalize = FALSE, 
    silent = FALSE, fun.kf = c("MARSSkfas", "MARSSkfss"))
}
\arguments{
  \item{object}{ An object of class \code{\link{marssMLE}}.}
  \item{method}{ Algorithm to use. Currently only "SS". }
  \item{normalize}{ TRUE/FALSE See details. }
  \item{silent}{ If TRUE, don't print inversion warnings. }
  \item{fun.kf}{ Can be ignored. This will change the Kalman filter/smoother function from the value in object$fun.kf if desired. }
}
\value{
A list with the following components  
  \item{model.residuals}{ The the observed one-step-ahead model residuals: data minus the model predictions conditioned on the data \eqn{t=1} to \eqn{t-1}. These are termed innovations. A n x T matrix. NAs will appear where the data are missing. }
  \item{state.residuals}{ The one-step-ahead state residuals \eqn{ \mathbf{x}_{t+1}^{t+1} - \mathbf{B}\mathbf{x}_{t}^t - \mathbf{u} }{ xtt(t+1) - B xtt(t) - u}
. Note, state residual at time \eqn{t} is the transition from time \eqn{t=t} to \eqn{t+1}. }
  \item{residuals}{ The residuals conditioned on the observed data up to time \eqn{t-1}. Returned as a (n+m) x T matrix with \code{model.residuals} in rows 1 to n and \code{state.residuals} in rows n+1 to n+m.  NAs will appear in rows 1 to n in the places where data are missing. }
  \item{var.residuals}{ The joint variance of the one-step-ahead residuals. Returned as a n+m x n+m x T matrix. }
  \item{std.residuals}{ The Cholesky standardized residuals as a n+m x T matrix. This is \code{residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals}. The model standardized residuals associated with the missing data are replaced with NA. }
  \item{mar.residuals}{ The marginal standardized residuals as a n+m x T matrix. This is \code{residuals} multiplied by the inverse of the diagonal matrix formed by the square-root of the diagonal of \code{var.residuals}.  The model marginal residuals associated with the missing data are replaced with NA. }
  \item{bchol.residuals}{ The Block Cholesky standardized residuals as a (n+m) x T matrix. This is \code{model.residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals[1:n,1:n,]} and \code{state.residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals[(n+1):(n+m),(n+1):(n+m),]}.}
  \item{E.obs.residuals}{ The expected value of the model residuals conditioned on the observed data \eqn{t=1} to \eqn{t-1}. Returned as a n x T matrix. Because all the data at time \eqn{t} are unobserved for the purpose of estimation (since conditioning is from \eqn{t=1} to \eqn{t-1}), this will be all 0s (unlike the case where we condition on the data from \eqn{t=1} to \eqn{T} or to \eqn{t}). This and \code{var.obs.residuals} are included for completeness since they are returned for \code{\link{MARSSresiduals.tT}()}, but they are not relevant for one-step-ahead residuals. See the discussion there. }
    \item{var.obs.residuals}{ For one-step-ahead residuals, this will be the same as the 1:n, 1:n upper diagonal block in \code{var.residuals} since none of the \eqn{t} data affect the residuals at time \eqn{t} (the model residuals are conditioned only on the data up to \eqn{t-1}). This is different for smoothation residuals which are conditioned on the data from \eqn{t=1} to \eqn{T}. This and \code{E.obs.residuals} are included for completeness since they are returned for \code{\link{MARSSresiduals.tT}()}, but they are not relevant for one-step-ahead residuals. See the discussion there. Note, also included as a code check. They are computed differently, but \code{var.obs.residuals} and \code{var.residuals} should always be the same.}
  \item{msg}{ Any warning messages. This will be printed unless \code{object$control$trace = -1} (suppress all error messages). }

}
\details{

This function returns the conditional expected value (mean) and variance of the one-step-ahead residuals.  'conditional' means in this context, conditioned on the observed data up to time \eqn{t-1} and a set of parameters.  

\strong{Model residuals}

\eqn{\mathbf{v}_t}{v_t} is the difference between the data and the predicted data at time \eqn{t} given \eqn{\mathbf{x}_t}{x_t}:
\deqn{ \mathbf{v}_t = \mathbf{y}_t - \mathbf{Z} \mathbf{x}_t - \mathbf{a} - \mathbf{D}\mathbf{d}_t}{ v(t) = y(t) - Z x(t) - a - D d(t)}
The observed model residuals \eqn{\hat{\mathbf{v}}_t}{hatv(t)} are the difference between the observed data and the predicted data at time \eqn{t} using the fitted model. \code{MARSSresiduals.tt1} fits the model using the data up to time \eqn{t-1}. So
\deqn{ \hat{\mathbf{v}}_t = \mathbf{y}_t - \mathbf{Z}\mathbf{x}_t^{t-1} - \mathbf{a}  - \mathbf{D}\mathbf{d}_t}{ hatv(t) = y(t) - Z xtt1(t) - a - D d(t)}
where \eqn{\mathbf{x}_t^{t-1}}{xtt1(t)} is the expected value of \eqn{\mathbf{X}_t}{X(t)} conditioned on the data from $t=1$ to \eqn{t-1} from the Kalman filter. \eqn{\mathbf{y}_t}{y(t)} are your data and missing values will appear as NA. 

\strong{State residuals}

\eqn{\mathbf{w}_{t+1}}{w(t+1)} are the difference between the state at time \eqn{t+1} and the expected value of the state at time \eqn{t+1} given the state at time \eqn{t}:
\deqn{ \mathbf{w}_{t+1} = \mathbf{x}_{t+1} - \mathbf{B} \mathbf{x}_{t} - \mathbf{u}  - \mathbf{C}\mathbf{c}_{t+1}}{ w(t+1) = x(t+1) - B x(t) - u - C c(t+1)}
The estimated state residuals \eqn{\hat{\mathbf{w}}_{t+1}}{hatw(t+1)} are the difference between estimate of \eqn{\mathbf{x}_{t+1}}{x(t+1)} minus the estimate using \eqn{\mathbf{x}_{t}}{x(t)}. 
\deqn{ \hat{\mathbf{w}}_{t+1} = \mathbf{x}_{t+1}^{t+1} - \mathbf{B}\mathbf{x}_{t}^t - \mathbf{u} - \mathbf{C}\mathbf{c}_{t+1}}{ hatw(t+1) = xtt(t+1) - B xtt(t) - u - C c(t+1)}
where \eqn{\mathbf{x}_{t+1}^{t+1}}{xtt(t+1)} is the Kalman filter estimate of the states at time \eqn{t+1} conditioned on the data up to time \eqn{t+1} and \eqn{\mathbf{x}_{t}^t}{xtt(t)} is the Kalman filter estimate of the states at time \eqn{t} conditioned on the data up to time \eqn{t}.
The estimated state residuals \eqn{\mathbf{w}_{t+1}}{w(t+1)} are returned in \code{state.residuals} and rows \eqn{n+1} to \eqn{n+m} of \code{residuals}. \code{state.residuals[,t]} is \eqn{\mathbf{w}_{t+1}}{w(t+1)} (notice time subscript difference). There are no NAs in the estimated state residuals (except for the last time step) as an estimate of the state exists whether or not there are associated data.

\code{res1} and \code{res2} in the code below will be the same.
\preformatted{dat <- t(harborSeal)[2:3,]
TT <- ncol(dat)
fit <- MARSS(dat)
B <- coef(fit, type="matrix")$B
U <- coef(fit, type="matrix")$U
xt <- MARSSkfss(fit)$xtt[,1:(TT-1)] # t 1 to TT-1
xtp1 <- MARSSkfss(fit)$xtt[,2:TT] # t 2 to TT
res1 <- xtp1 - B \%*\% xt - U \%*\% matrix(1,1,TT-1)
res2 <- MARSSresiduals(fit, type="tt1")$state.residuals
}

\strong{Joint residual variance}

In a state-space model, \eqn{\mathbf{X}}{X} and \eqn{\mathbf{Y}}{Y} are stochastic, and the model and state residuals are random variables \eqn{\hat{\mathbf{V}}_t}{hatV(t)} and \eqn{\hat{\mathbf{W}}_{t+1}}{hatW(t+1)}. The joint distribution of \eqn{\hat{\mathbf{V}}_{t}, \hat{\mathbf{W}}_{t+1}}{hatV(t), hatW(t+1)} is the distribution across all the different possible data sets that our MARSS equations with parameters \eqn{\Theta}{Theta} might generate. Denote the matrix of \eqn{\hat{\mathbf{V}}_{t}, \hat{\mathbf{W}}_{t+1}}{hatV(t), hatW(t+1)}, as \eqn{\widehat{\mathcal{E}}_{t}}{Epsilon(t)}. That distribution has an expected value (mean) and variance:
\deqn{ \textrm{E}[\widehat{\mathcal{E}}_t] = 0; \textrm{var}[\widehat{\mathcal{E}}_t] = \hat{\Sigma}_t }{ E[Epsilon(t)] = 0; var[Epsilon(t)] = hatSigma(t)}
Our observed residuals \code{residuals} are one sample from this distribution.
To standardize the observed residuals, we will use \eqn{ \hat{\Sigma}_t }{ hatSigma(t) }. \eqn{ \hat{\Sigma}_t }{ hatSigma(t) } is returned in \code{var.residuals}. Rows/columns 1 to \eqn{n} are the conditional variances of the model residuals and rows/columns \eqn{n+1} to \eqn{n+m} are the conditional variances of the state residuals. The off-diagonal blocks are the covariances between the two types of residuals. For one-step-ahead residuals (unlike smoothation residuals \link{MARSSresiduals.tT}), the covariance is zero.

\code{var.residuals} returned by this function is the conditional variance of the residuals conditioned on the data up to \eqn{t-1} and the parameter set \eqn{\Theta}{Theta}.  The conditional variance for the model residuals is 
\deqn{ \hat{\Sigma}_t = \mathbf{R}+\mathbf{Z}_t \mathbf{V}_t^{t-1} \mathbf{Z}_t^\top }{hatSigma(t) = R + Z Vtt1 t(Z)}
where \eqn{\mathbf{V}_t^{t-1}}{Vtt1} is the variance of \eqn{\mathbf{X}_t}{X(t)} conditioned on the data up to time \eqn{t-1}. This is returned by \code{\link{MARSSkf}} in \code{Vtt1}. The innovations variance is also returned in \code{Sigma} from \code{\link{MARSSkf}} and are used in the innovations form of the likelihood calculation.

\strong{Standardized residuals}

\code{std.residuals} are Cholesky standardized residuals. These are the residuals multiplied by the inverse of the lower triangle of the Cholesky decomposition of the variance matrix of the residuals: 
\deqn{ \hat{\Sigma}_t^{-1/2} \hat{\mathbf{v}}_t}{ hatSigma(t)^{-1/2} hatv(t). }
These residuals are uncorrelated unlike marginal residuals. 

The interpretation of the Cholesky standardized residuals is not straight-forward when the \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R} variance-covariance matrices are non-diagonal.  The residuals which were generated by a non-diagonal variance-covariance matrices are transformed into orthogonal residuals in MVN(0,I) space.  For example, if v is 2x2 correlated errors with variance-covariance matrix R. The transformed residuals (from this function) for the i-th row of v is a combination of the row 1 effect and the row 1 effect plus the row 2 effect.  So in this case, row 2 of the transformed residuals would not be regarded as solely the row 2 residual but rather how different row 2 is from row 1, relative to expected.  If the errors are highly correlated, then the Cholesky standardized residuals can look rather non-intuitive.

\code{mar.residuals} are the marginal standardized residuals. These are the residuals multiplied by the inverse of the diagonal matrix formed from the square-root of the diagonal of the variance matrix of the residuals: 
\deqn{ \textrm{dg}(\hat{\Sigma}_t)^{-1/2} \hat{\mathbf{v}}_t}{ dg(hatSigma(t))^{-1/2} hatv(t)}, where 'dg(A)' is the square matrix formed from the diagonal of A, aka \code{diag(diag(A))}. These residuals will be correlated if the variance matrix is non-diagonal. 

The Block Cholesky standardized residuals are like the Cholesky standardized residuals except that the full variance-covariance matrix is not used, only the variance-covariance matrix for the model or state residuals (respectively) is used for standardization. For the one-step-ahead case, the model and state residuals are independent (unlike in the smoothations case) thus the Cholesky and Block Cholesky standardized residuals will be identical (unlike in the smoothations case). 

\strong{Normalized residuals}

If \code{normalize=FALSE}, the unconditional variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)} are \eqn{\mathbf{R}}{R} and \eqn{\mathbf{Q}}{Q} and the model is assumed to be written as
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{v}_t}{ y(t) = Z x(t) + a + v(t)}
\deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{w}_t}{ x(t) = B x(t-1) + u + w(t)}
If normalize=TRUE, the model is assumed to be written
\deqn{\mathbf{y}_t = \mathbf{Z} \mathbf{x}_t + \mathbf{a} + \mathbf{H}\mathbf{v}_t}{ y(t) = Z x(t) + a + Hv(t)}
\deqn{\mathbf{x}_t = \mathbf{B} \mathbf{x}_{t-1} + \mathbf{u} + \mathbf{G}\mathbf{w}_t}{ x(t) = B x(t-1) + u + Gw(t)}
with the variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)} equal to \eqn{\mathbf{I}}{I} (identity).

\code{MARSSresiduals} returns the residuals defined as in the first equations. To get the residuals defined as Harvey et al. (1998) define them (second equations), then use \code{normalize=TRUE}.  In that case the unconditional variance of residuals will be \eqn{\mathbf{I}}{I} instead of \eqn{\mathbf{Q}}{Q} and \eqn{\mathbf{R}}{R}.  Note, that the normalized residuals are not the same as the standardized residuals.  In former, the unconditional residuals have a variance of \eqn{\mathbf{I}}{I} while in the latter it is the conditional residuals that have a variance of \eqn{\mathbf{I}}{I}.

}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
} 
\seealso{ \code{\link{MARSSresiduals.tT}()}, \code{\link{MARSSresiduals.tt}()},  \code{\link{fitted.marssMLE}()}, \code{\link{plot.marssMLE}()} }
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  fit <- MARSS(dat)
  
  MARSSresiduals(fit, type="tt1")$std.residuals
  residuals(fit, type="tt1")
}
\references{
R. H. Shumway and D. S. Stoffer (2006).  Section on the calculation of the likelihood of state-space models in Time series analysis and its applications.  Springer-Verlag, New York.

Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045. 
}

