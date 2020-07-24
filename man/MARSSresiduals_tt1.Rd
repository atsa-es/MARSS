\name{MARSSresiduals.tt1}
\alias{MARSSresiduals.tt1}

\title{ MARSS One-Step-Ahead Residuals }
\description{
  Calculates the standardized (or auxiliary) one-step-ahead residuals, aka the innovations residuals and their variance.  One-step-ahead residuals are only for the observatons. Not exported. Access this function with \code{MARSSresiduals(object, type="innovations")}.
}
\usage{
MARSSresiduals.tt1(object, method=c("SS"), normalize=FALSE, silent=FALSE)
}
\arguments{
  \item{object}{ An object of class \code{\link{marssMLE}}.}
  \item{method}{ Algorithm to use. Currently only "SS". }
  \item{normalize}{ TRUE/FALSE See details. }
  \item{silent}{ If TRUE, don't print inversion warnings. }
}
\value{
A list with the following components  
  \item{model.residuals}{ The the observed one-step-ahead model residuals: data minus the model predictions conditioned on the data 1 to t-1. These are termed innovations. A n x T matrix. NAs will appear where the data are missing. }
  \item{state.residuals}{ NULL There are no one-step-ahead residuals for the states. }
  \item{residuals}{ The residuals. Since there are only model residuals this is the same as \code{model.residuals}. }
  \item{var.residuals}{ The joint variance of the one-step-ahead model residuals conditioned on observed data from 1 to t-1. Returned as a n x n x T matrix. }
  \item{std.residuals}{ The Cholesky standardized residuals as a n x T matrix. This is \code{residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals}. The model standardized residuals associated with the missing data are replaced with NA. }
  \item{mar.residuals}{ The marginal standardized residuals as a n x T matrix. This is \code{residuals} multiplied by the inverse of the diagonal matrix formed by the square-root of the diagonal of \code{var.residuals}.  The model marginal residuals associated with the missing data are replaced with NA. }
    \item{E.obs.residuals}{ The expected value of the model residuals conditioned on the observed data 1 to t-1. Returned as a n x T matrix. Because all the data at time t are unobserved for the purpose of estimation (since conditioning is from 1 to t-1), this will be all 0s (unlike the case where we condition on the data from 1 to T). This and \code{var.obs.residuals} are included for completeness since they are returned for \code{\link{MARSSresiduals.tT}}, but they are not relevant for one-step-ahead residuals. See the discussion there. }
    \item{var.obs.residuals}{ For one-step-ahead residuals, this will be the same as \code{var.residuals} since none of the t data affect the residual at time t, since the residuals are conditioned only on the data up to t-1. This is different for residuals conditioned on the data from 1 to T. This and \code{E.obs.residuals} are included for completeness since they are returned for \link{\code{MARSSresiduals.tT}}, but they are not relevant for one-step-ahead residuals. See the discussion there. Note, also included as a code check. They are computed differently, but \code{var.obs.residuals} and \code{var.residuals} should always be the same.}
  \item{msg}{ Any warning messages. This will be printed unless Object$control$trace = -1 (suppress all error messages). }

}
\details{

This function returns the conditional expected value (mean) and variance of the model one-step-ahead residuals.  'conditional' means in this context, conditioned on the observed data up to time \eqn{t-1} and a set of parameters.  

\strong{Model residuals}

\eqn{\mathbf{v}_t}{v_t} is the difference between the data and the predicted data at time \eqn{t} given \eqn{\mathbf{x}_t}{x_t}:
\deqn{ \mathbf{v}_t = \mathbf{y}_t - \mathbf{Z} \mathbf{x}_t - \mathbf{a} }{ v_t = y_t - Z x_t - a}
The observed model residuals \eqn{\hat{\mathbf{v}}_t}{hatv_t} are the difference between the observed data and the predicted data at time \eqn{t} using the fitted model. \code{MARSSresiduals.tt1} fits the model using the data up to time \eqn{t-1}. So
\deqn{ \hat{\mathbf{v}}_t = \mathbf{y}_t - \mathbf{Z}\tilde{\mathbf{x}}_t^{t-1} - \mathbf{a} }{ hatv_t = y_t - Z xtt1 - a}
where \eqn{\tilde{\mathbf{x}}_t^{t-1}}{xtt1} is the expected value of \eqn{\mathbf{X}_t}{X_t} conditioned on the data from 1 to \eqn{t-1} from the Kalman filter. \eqn{\mathbf{y}_t}{y_t} are your data and missing values will appear as NA. These will be returned in \code{residuals}.

\code{var.residuals} returned by the function is the conditional variance of the residuals conditioned on the data up to \eqn{t-1} and the parameter set \eqn{\Theta}.  The conditional variance is 
\deqn{ \hat{\Sigma}_t = \mathbf{R}+\mathbf{Z}_t \mathbf{V}_t^{t-1} \mathbf{Z}_t^\top }{hatSigma_t = R + Z_t Vtt1 t(Z_t)}
where \eqn{\mathbf{V}_t^{t-1}}{Vtt1} is the variance of \eqn{\mathbf{X}_t}{X_t} conditioned on the data up to time \eqn{t-1}. This is returned by \code{\link{MARSSkf}} in \code{Vtt1}.

\strong{Standardized residuals}

\code{std.residuals} are Cholesky standardized residuals. These are the residuals muliplied by the inverse of the lower triangle of the Cholesky decomposition of the variance matrix of the residuals: 
\deqn{ \hat{\Sigma}_t^{-1/2} \hat{\mathbf{v}}_t}{ hatSigma_t^{-1/2} hatv_t. }
These residuals are uncorrelated unlike marginal residuals. 

The interpretation of the Cholesky standardized residuals is not straight-forward when the Q and R variance-covariance matrices are non-diagonal.  The residuals which were generated by a non-diagonal variance-covariance matrices are transformed into orthogonal residuals in MVN(0,I) space.  For example, if v is 2x2 correlated errors with variance-covariance matrix R. The transformed residuals (from this function) for the i-th row of v is a combination of the row 1 effect and the row 1 effect plus the row 2 effect.  So in this case, row 2 of the transformed residuals would not be regarded as solely the row 2 residual but rather how different row 2 is from row 1, relative to expected.  If the errors are highly correlated, then the Cholesky standardized residuals can look rather non-intuitive.

\code{mar.residuals} are the marginal standardized residuals. These are the residuals muliplied by the inverse of the diagonal matrix formed from the square-root of the diagonal of the variance matrix of the residuals: 
\deqn{ \textrm{dg}(\hat{\Sigma}_t)^{-1/2} \hat{v}_t}{ dg(hatSigma_t)^{-1/2} hatv_t}, where 'dg(A)' is the square matrix formed from the diagonal of A, aka \code{diag(diag(A))}. These residuals will be correlated if the variance matrix is non-diagonal. 

\strong{Normalized residuals}

If \code{normalize=FALSE}, the unconditional variance of \eqn{V_t} and \eqn{W_t} are R and Q and the model is assumed to be written as
\deqn{ y_t = Z x_t + a + v_t}
\deqn{ x_t = B x_{t-1} + u + w_t}
If normalize=TRUE, the model is assumed to be written
\deqn{ y_t = Z x_t + a + Hv_t}
\deqn{ x_t = B x_{t-1} + u + Gw_t}
with the variance of \eqn{V_t} and \eqn{W_t} equal to I (identity).

\code{MARSSresiduals} returns the residuals defined as in the first equations. To get the residuals defined as Harvey et al. (1998) define them (second equations), then use \code{normalize=TRUE}.  In that case the unconditional variance of residuals will be I instead of R and Q.  Note, that the `normalized' residuals are not the same as the `standardized' residuals.  In former, the unconditional residuals have a variance of I while in the latter it is the conditional residuals that have a variance of I.

}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSSresiduals.tT}}, \code{\link{fitted.marssMLE}}, \code{\link{plot.marssMLE}} }
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  MLEobj <- MARSS(dat)
  
  MARSSresiduals(MLEobj, type="tt1")$std.residuals
  residuals(MLEobj, type="innovations")
}
\references{
R. H. Shumway and D. S. Stoffer (2006).  Section on the calculation of the likelihood of state-space models in Time series analysis and its applications.  Springer-Verlag, New York.

Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045. 
}

