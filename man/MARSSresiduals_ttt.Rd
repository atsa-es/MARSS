\name{MARSSresiduals.tt}
\alias{MARSSresiduals.tt}

\title{ MARSS Contemporaneous Residuals }
\description{
  Calculates the standardized (or auxiliary) contemporaneous residuals, aka the residuals and their variance conditioned on the data up to time \eqn{t}.  Contemporaneous residuals are only for the observations. Not exported. Access this function with \code{MARSSresiduals(object, type="tt")}.
}
\usage{
MARSSresiduals.tt(object, method=c("SS"), normalize=FALSE, 
    silent=FALSE, fun.kf=c("MARSSkfas", "MARSSkfss"))
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
  \item{model.residuals}{ The the observed contemporaneous model residuals: data minus the model predictions conditioned on the data 1 to t. A n x T matrix. NAs will appear where the data are missing. }
  \item{state.residuals}{ All NA. There are no contemporaneous residuals for the states. }
  \item{residuals}{ The residuals. \code{model.residuals} are in rows 1:n and \code{state.residuals} are in rows n+1:n+m. }
  \item{var.residuals}{ The joint variance of the residuals conditioned on observed data from 1 to t-. This only has values in the 1:n,1:n upper block for the model residuals. }
  \item{std.residuals}{ The Cholesky standardized residuals as a n+m x T matrix. This is \code{residuals} multiplied by the inverse of the lower triangle of the Cholesky decomposition of \code{var.residuals}. The model standardized residuals associated with the missing data are replaced with NA. Note because the contemporaneous state residuals do not exist, rows n+1:n+m are all NA. }
  \item{mar.residuals}{ The marginal standardized residuals as a n+m x T matrix. This is \code{residuals} multiplied by the inverse of the diagonal matrix formed by the square-root of the diagonal of \code{var.residuals}.  The model marginal residuals associated with the missing data are replaced with NA. }
  \item{bchol.residuals}{ Because state residuals do not exist, this will be equivalent to the Cholesky standardized residuals, \code{std.residuals}. }
  \item{E.obs.residuals}{ The expected value of the model residuals conditioned on the observed data 1 to t. Returned as a n x T matrix. }
    \item{var.obs.residuals}{ The variance of the model residuals conditioned on the observed data. Returned as a n x n x T matrix.  For observed data, this will be 0. See  \code{\link{MARSSresiduals.tT}} for a discussion of these residuals and where they might be used. }
  \item{msg}{ Any warning messages. This will be printed unless Object$control$trace = -1 (suppress all error messages). }

}
\details{

This function returns the conditional expected value (mean) and variance of the model contemporaneous residuals.  'conditional' means in this context, conditioned on the observed data up to time \eqn{t} and a set of parameters.  

\strong{Model residuals}

\eqn{\mathbf{v}_t}{v(t)} is the difference between the data and the predicted data at time \eqn{t} given \eqn{\mathbf{x}_t}{x(t)}:
\deqn{ \mathbf{v}_t = \mathbf{y}_t - \mathbf{Z} \mathbf{x}_t - \mathbf{a} - \mathbf{d}\mathbf{d}_{t}}{ v(t) = y(t) - Z x(t) - a - D d(t)}
The observed model residuals \eqn{\hat{\mathbf{v}}_t}{hatv(t)} are the difference between the observed data and the predicted data at time \eqn{t} using the fitted model. \code{MARSSresiduals.tt} fits the model using the data up to time \eqn{t}. So
\deqn{ \hat{\mathbf{v}}_t = \mathbf{y}_t - \mathbf{Z}\mathbf{x}_t^{t} - \mathbf{a} - \mathbf{D}\mathbf{d}_{t}}{ hatv(t) = y(t) - Z xtt - a - D d(t)}
where \eqn{\mathbf{x}_t^{t}}{xtt(t)} is the expected value of \eqn{\mathbf{X}_t}{X(t)} conditioned on the data from 1 to \eqn{t} from the Kalman filter. \eqn{\mathbf{y}_t}{y(t)} are your data and missing values will appear as NA. These will be returned in \code{residuals}.

\code{var.residuals} returned by the function is the conditional variance of the residuals conditioned on the data up to \eqn{t} and the parameter set \eqn{\Theta}{Theta}.  The conditional variance is 
\deqn{ \hat{\Sigma}_t = \mathbf{R}+\mathbf{Z} \mathbf{V}_t^{t} \mathbf{Z}^\top }{hatSigma(t) = R + Z Vtt t(Z)}
where \eqn{\mathbf{V}_t^{t}}{Vtt} is the variance of \eqn{\mathbf{X}_t}{X(t)} conditioned on the data up to time \eqn{t}. This is returned by \code{\link{MARSSkfss}} in \code{Vtt}.

\strong{Standardized residuals}

\code{std.residuals} are Cholesky standardized residuals. These are the residuals multiplied by the inverse of the lower triangle of the Cholesky decomposition of the variance matrix of the residuals: 
\deqn{ \hat{\Sigma}_t^{-1/2} \hat{\mathbf{v}}_t}{ hatSigma(t)^{-1/2} hatv(t) }.
These residuals are uncorrelated unlike marginal residuals. 

The interpretation of the Cholesky standardized residuals is not straight-forward when the Q and R variance-covariance matrices are non-diagonal.  The residuals which were generated by a non-diagonal variance-covariance matrices are transformed into orthogonal residuals in MVN(0,I) space.  For example, if v is 2x2 correlated errors with variance-covariance matrix R. The transformed residuals (from this function) for the i-th row of v is a combination of the row 1 effect and the row 1 effect plus the row 2 effect.  So in this case, row 2 of the transformed residuals would not be regarded as solely the row 2 residual but rather how different row 2 is from row 1, relative to expected.  If the errors are highly correlated, then the Cholesky standardized residuals can look rather non-intuitive.

\code{mar.residuals} are the marginal standardized residuals. These are the residuals multiplied by the inverse of the diagonal matrix formed from the square-root of the diagonal of the variance matrix of the residuals: 
\deqn{ \textrm{dg}(\hat{\Sigma}_t)^{-1/2} \hat{v}_t}{ dg(hatSigma(t))^{-1/2} hatv(t)}, where 'dg(A)' is the square matrix formed from the diagonal of A, aka \code{diag(diag(A))}. These residuals will be correlated if the variance matrix is non-diagonal. 

\strong{Normalized residuals}

If \code{normalize=FALSE}, the unconditional variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)} (the model and state residuals respectively) are \eqn{\mathbf{R}}{R} and \eqn{\mathbf{Q}}{Q} and the model is assumed to be written as
\deqn{ y_t = Z x_t + a + v_t}
\deqn{ x_t = B x_{t-1} + u + w_t}
If normalize=TRUE, the model is assumed to be written
\deqn{ y_t = Z x_t + a + Hv_t}
\deqn{ x_t = B x_{t-1} + u + Gw_t}
with the variance of \eqn{\mathbf{V}_t}{V(t)} and \eqn{\mathbf{W}_t}{W(t)} equal to I (identity).

\code{MARSSresiduals} returns the residuals defined as in the first equations. To get normalized residuals (second equation) as used in Harvey et al. (1998), then use \code{normalize=TRUE}.  In that case the unconditional variance of residuals will be I instead of \eqn{\mathbf{R}}{R} and \eqn{\mathbf{Q}}{Q}.  Note, that the `normalized' residuals are not the same as the `standardized' residuals.  In former, the unconditional residuals have a variance of I while in the latter it is the conditional residuals that have a variance of I.

}

\author{ 
  Eli Holmes, NOAA, Seattle, USA.  

  eli(dot)holmes(at)noaa(dot)gov
} 
\seealso{ \code{\link{MARSSresiduals.tT}}, \code{\link{MARSSresiduals.tt1}}, \code{\link{fitted.marssMLE}}, \code{\link{plot.marssMLE}} }
\examples{
  dat <- t(harborSeal)
  dat <- dat[c(2,11),]
  fit <- MARSS(dat)
  
  # Returns a matrix
  MARSSresiduals(fit, type="tt")$std.residuals
  # Returns a data frame in long form
  residuals(fit, type="tt")
}
\references{
Holmes, E. E. 2014. Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045. 
}

