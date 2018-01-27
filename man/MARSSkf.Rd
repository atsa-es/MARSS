\name{MARSSkf}
\alias{MARSSkf}
\alias{MARSSkfas}
\alias{MARSSkfss}
\title{ Kalman Filtering and Smoothing for Time-varying MARSS models}
\description{
  Provides Kalman filter and smoother output for MARSS models with (or without) time-varying parameters. This is a base function in the \code{\link{MARSS-package}}.  \code{MARSSkf} is a small helper function to select which Kalman filter/smoother function to use based on which function was requested (in \code{MLEobj$fun.kf}).  The default function is \code{MARSSkfas}.
}
\usage{
MARSSkf( MLEobj, only.logLik=FALSE, return.lag.one=TRUE, return.kfas.model=FALSE )
MARSSkfss( MLEobj )
MARSSkfas( MLEobj, only.logLik=FALSE, return.lag.one=TRUE, return.kfas.model=FALSE )
}
\arguments{
  \item{ MLEobj }{ A \code{\link{marssMLE}} object with the \code{par} element of estimated parameters, \code{marss} element with the model description (in marss form) and data, and \code{control} element for the fitting algorithm specificitions.  \code{control$debugkf} specifies that detailed error reporting will be returned (only used by \code{MARSSkf}).  \code{model$diffuse=TRUE} specifies that a diffuse prior be used (only used by \code{MARSSkfas}). See KFAS documentation. When the diffuse prior is set, V0 should be non-zero since the diffuse prior variance is V0*kappa, where kappa goes to infinity.}
  \item{ only.logLik }{ Used by \code{MARSSkfas}.  If set, only the log-likelihood is returned using the KFAS function \code{\link[KFAS]{logLik.SSModel}}.  This is much faster if only the log-likelihood is needed. }
  \item{ return.lag.one }{ Used by \code{MARSSkfas}.  If set to FALSE, the smoothed lag-one covariance  values are not returned (Vtt1T is set to NULL).  This speeds up \code{MARSSkfas} because to return the smoothed lag-one covariance a stacked MARSS model is used with twice the number of state vectors---thus the state matrices are larger and take more time to work with. }
  \item{ return.kfas.model }{ Used by \code{MARSSkfas}.  If set to TRUE, it returns the MARSS model in KFAS model form (class \code{\link[KFAS]{SSModel}}).  This is useful if you want to use other KFAS functions or write your own functions to work with \code{\link{optim}} to do optimization.  This can speed things up since there is a bit of code overhead in \code{\link{MARSSoptim}} associated with the \code{\link{marssMODEL}} model specification needed for the constrained EM algorithm (but not strictly needed for \code{\link{optim}}; useful but not required.). }
}
\details{
For state-space models, the Kalman filter and smoother provide optimal (minimum mean square error) estimates of the hidden states. The Kalman filter is a forward recursive algorithm which computes estimates of the states x(t) conditioned on the data up to time t (xtt). The Kalman smoother is a backward recursive algorithm which starts at time T and works backwards to t = 1 to provide estimates of the states conditioned on all data (xtT).    The data may contain missing values (NAs).  All parameters may be time varying.

The expected value of the initial state, x0, is an estimated parameter (or treated as a prior).  This E(initial state) can be treated in two different ways.  One can treat it as x00, meaning E(x at t=0 | y at t=0), and then compute x10, meaning E(x at t=1 | y at t=0), from x00.    Or one can simply treat the initial state as x10, meaning E(x at t=1 | y at t=0).  The approaches lead to the same parameter estimates, but the likelihood is written slightly differently in each case and you need your likelihood calculation to correspond to how the initial state is treated in your model (either x00 or x10).  The EM algorithm in the MARSS package (\code{\link{MARSSkem}}) provides both Shumway and Stoffer's derivation that uses tinitx=0 and Ghahramani et al algorithm which uses tinitx=1.  The \code{MLEobj$model$tinitx} argument specifies whether the initial states (specified with x0 and V0) is at t=0 (\code{tinitx=0}) or t=1 (\code{tinitx=1}).

\code{MARSSkfss()} is a native R implementation based on the traditional Kalman filter and smoother equation as shown in Shumway and Stoffer (sec 6.2, 2006).  The equations have been altered slightly to the initial state distribution to be to be specified at t=0 or t=1 (data starts at t=1) per per Ghahramani and Hinton (1996).  In addition, the filter and smoother equations have been altered to allow partially deterministic models (some or all elements of the Q diagonals equal to 0), partially perfect observation models (some or all elements of the R diagonal equal to 0) and fixed (albeit unknown) initial states (some or all elements of the V0 diagonal equal to 0) (per Holmes 2012).  The code includes numerous checks to alert the user if matrices are becoming ill-conditioned and the algorithm unstable.  

\code{MARSSkfas()} uses the (Fortran-based) Kalman filter and smoother function (\code{\link[KFAS]{KFS}}) in the KFAS package (Helske 2012) based on the algorithms of Koopman and Durbin (2000, 2001, 2003).  The Koopman and Durbin algorithm is faster and more stable since it avoids matrix inverses.  Exact diffuse priors are also allowed in the KFAS Kalman filter function.  The standard output from the KFAS functions do not include the lag-one covariance smoother needed for the EM algorithm.  \code{MARSSkfas} computes the smoothed lag-one covariance  using the Kalman filter applied to a stacked MARSS model as described on page 321 in Shumway and Stoffer (2000). Also the KFAS model specification only has the initial state at t=1 (as x(1) conditioned on y(0), which is missing).  When the initial state is specified at t=0 (as x(0) conditioned on y(0), which is missing), \code{MARSSkfas} computes the required E(x(1)|y(0)) and var(x(1)|y(0)) using the Kalman filter equations per Ghahramani and Hinton (1996).  
 
The likelihood returned for both functions is the exact likelihood when there are missing values rather than the approximate likelihood sometimes presented in texts for the missing values case.  The functions return the same filter, smoother and log-likelihood values.  The differences are that \code{MARSSkfas} is faster (and more stable) but \code{MARSSkf} has many internal checks and error messages which can help debug numerical problems (but slow things down).  Also  \code{MARSSkf} returns some output specific to the traditional filter algorithm (J and Kt).  
}
\value{
A list with the following components (m is the number of state processes). "V" elements are called "P" in Shumway and Stoffer (S&S eqn 6.17 with s=T).  The output is referenced against equations in Shumway and Stoffer (2006); the Kalman filter and smoother implemented in MARSS is for a more general MARSS model than that shown in S&S but the output has the same meaning.  In the expectations below, the parameters are left off, so E[x | y] is really E[x | theta, y] where theta is the parameter list.
  \item{xtT}{ State first moment conditioned on y(1:T): E[x(t) | y(1:T)] (m x T matrix). Kalman smoother output.}
  \item{VtT}{ State variance conditioned on y(1:T): E[(x(t)-xtT(t))(x(t)-xtT(t))'| | y(1:T)] (m x m x T array). Kalman smoother output. P_t^T in S&S (S&S eqn 6.18 with s=T, t1=t2=t). State second moment E[x(t)x(t)'| y(1:T)] = VtT(t)+xtT(t)xtT(t)'}
  \item{Vtt1T}{ State lag-one covariance E[(x(t)-xtT(t))(x(t-1)-xtT(t-1))' | y(1:T)] (m x m x T). Kalman smoother output. P_{t,t-1}^T in S&S (S&S eqn 6.18 with s=T, t1=t, t2=t-1). State lag-one second moments E[x(t)x(t-1)'| y(1:T)] = Vtt1T(t)+xtT(t)xtT(t-1)'. }
  \item{x0T}{ Initial state estimate E[x(i) | y(1:T)] (m x 1). If control$kf.x0="x00", i=0; if ="x10", i=1. Kalman smoother output.}
  \item{V0T}{ Estimate of initial state covariance matrix E[x(i)x(i)' | y(1:T)] (m x m). If model$tinitx=0, i=0; if =1, i=1.  Kalman smoother output. P_0^T in S&S.}
  \item{J}{ (m x m x T) Kalman smoother output.  Only for \code{MARSSkfss}. (ref S&S eqn 6.49)}
  \item{J0}{ J at init time (t=0 or t=1) (m x m x T). Kalman smoother output. Only for \code{MARSSkfss}.}
  \item{xtt}{ State first moment conditioned on y(1:t): E[x(t) | y(1:t)] (m x T).  Kalman filter output. (S&S eqn 6.17 with s=t)}
  \item{xtt1}{ State first moment conditioned on y(1:t-1): E[x(t) | y(1:t-1)] (m x T).  Kalman filter output. (S&S eqn 6.17 with s=t-1) }
  \item{Vtt}{ State variance conditioned on y(1:t): E[(x(t)-xtt(t))(x(t)-xtt(t))'| | y(1:t)] (m x m x T array). Kalman filter output. P_t^t in S&S (S&S eqn 6.18 with s=t, t1=t2=t). State second moment E[x(t)x(t)'| y(1:t)] = Vtt(t)+xtt(t)xtt(t)' }
  \item{Vtt1}{ State variance conditioned on y(1:t-1): E[(x(t)-xtt1(t))(x(t)-xtt1(t))'| | y(1:t-1)] (m x m x T array). Kalman filter output. P_t^t in S&S (S&S eqn 6.18 with s=t-1, t1=t2=t). State second moment E[x(t)x(t)'| y(1:t-1)] = Vtt1(t)+xtt1(t)xtt1(t)'}
  \item{Kt}{ Kalman gain (m x m x T). Kalman filter output (ref S&S eqn 6.23). Only for \code{MARSSkfss}.}
  \item{Innov}{ Innovations y(t) - E[y(t) | y(1:t-1)] (n x T). Kalman filter output. Only returned with \code{MARSSkfss}. (ref page S&S 339). }
  \item{Sigma}{ Innovations covariance matrix. Kalman filter output. Only returned with \code{MARSSkfss}. (ref S&S eqn 6.61)}
  \item{logLik}{ Log-likelihood logL(y(1:T) | theta) (ref S&S eqn 6.62)}
  \item{kfas.model}{ The model in KFAS model form (class \code{\link[KFAS]{SSModel}}). Only for \code{MARSSkfas}.}
  \item{errors}{ Any error messages. }
}
\references{
A. C. Harvey (1989).  Chapter 5, Forecasting, structural time series models and the Kalman filter.  Cambridge University Press.

R. H. Shumway and D. S. Stoffer (2006).  Time series analysis and its applications: with R examples.  Second Edition. Springer-Verlag, New York.

Ghahramani, Z. and Hinton, G.E. (1996) Parameter estimation for linear dynamical systems.  University of Toronto Technical Report CRG-TR-96-2.

Holmes, E. E. (2012).  Derivation of the EM algorithm for constrained and unconstrained multivariate autoregressive
state-space (MARSS) models.  Technical Report. arXiv:1302.3919 [stat.ME] \code{RShowDoc("EMDerivation",package="MARSS")} to open a copy.

Jouni Helske (2012). KFAS: Kalman filter and smoother for exponential family state space models. R package version
  0.9.11. \code{http://CRAN.R-project.org/package=KFAS}
  
Koopman, S.J. and Durbin J. (2000). Fast filtering and smoothing for non-stationary time series models, Journal of American Statistical Assosiation, 92, 1630-38.

Koopman, S.J. and Durbin J. (2001). Time series analysis by state space methods. Oxford: Oxford University Press.

Koopman, S.J. and Durbin J. (2003). Filtering and smoothing of state vector for diffuse state space models, Journal of Time Series Analysis, Vol. 24, No. 1.

The user guide:  Holmes, E. E., E. J. Ward, and M. D. Scheuerell (2012) Analysis of multivariate time-series using the MARSS package. NOAA Fisheries, Northwest Fisheries Science Center, 2725 Montlake Blvd E., Seattle, WA 98112   Type \code{RShowDoc("UserGuide",package="MARSS")} to open a copy.     

}
\author{ 
  Eli Holmes, NOAA, Seattle, USA.  
  eli(dot)holmes(at)noaa(dot)gov
}
\seealso{ 
\code{\link{MARSS}} \code{\link{marssMODEL}} \code{\link{MARSSkem}}
}

\examples{
  dat = t(harborSeal)
  dat = dat[2:nrow(dat),]
  #you can use MARSS to construct a MLEobj
  #MARSS calls MARSSinits to construct default initial values
  MLEobj = MARSS(dat, fit=FALSE)
  #MARSSkf needs a marss MLE object with the par element set
  MLEobj$par=MLEobj$start
  #Compute the kf output at the params used for the inits 
  kfList = MARSSkfas( MLEobj )
}

