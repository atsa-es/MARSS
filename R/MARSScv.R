#' MARSScv is a wrapper for MARSS that re-fits the model with cross validated data.
#'
#' @param y A n x T matrix of n time series over T time steps. Only y
#' is required for the function. A ts object (univariate or multivariate)
#' can be used and this will be converted to a matrix with time in the
#' columns.
#' @param model Model specification using a list of parameter matrix text shortcuts or matrices.
#' See Details and \code{\link{MARSS.marxss}} for the default form.
#' Or better yet open the Quick Start Guide \code{RShowDoc("Quick_Start",package="MARSS")}
#' @param inits A list with the same form as the list outputted by \code{coef(fit)}
#' that specifies initial values for the parameters. See also \code{\link{MARSS.marxss}}.
#' @param method Estimation method. MARSS provides an EM algorithm (\code{method="kem"})
#' (see \code{\link{MARSSkem}}) and the BFGS algorithm (\code{method="BFGS"})
#' (see \code{\link{MARSSoptim}})
#' @param form The equation form used in the \code{MARSS()} call.  The default is "marxss".
#' See \code{\link{MARSS.marxss}} or \code{\link{MARSS.dfa}}
#' @param silent Setting to TRUE(1) suppresses printing of full error messages, warnings,
#' progress bars and convergence information. Setting to FALSE(0) produces
#' error output. Setting silent=2 will produce more verbose error messages and
#' progress information.
#' @param control Estimation options for the maximization algorithm. The typically used
#' control options for method="kem" are below but see  \code{\link{marssMLE}} for the full
#' list of control options.  Note many of these are not allowed if method="BFGS";
#' see  \code{\link{MARSSoptim}} for the allowed control options for this method.
#' @param fun.kf What Kalman filter function to use.  MARSS has two:
#' \code{\link{MARSSkfas}()} which is based on the Kalman filter in the
#' \href{https://cran.r-project.org/package=KFAS}{KFAS} package based on
#' Koopman and Durbin and \code{\link{MARSSkfss}()} which is a native
#' R implementation of the Kalman filter and smoother in Shumway and
#' Stoffer. The KFAS filter is much faster.  \code{\link{MARSSkfas}()}
#' modifies the input and output in order to output the lag-one covariance
#' smoother needed for the EM algorithm (per page 321 in Shumway and Stoffer (2000).
#' @param fold_ids A n x T matrix of integers, with values assigned by the user to folds.
#' If not included, data are randomly assigned to one of 10 folds
#' @param future_cv Whether or not to use future cross validation (defaults to FALSE), where data up to
#' time T-1 are used to predict data at time T. Data are held out by time slices, and
#' the `fold_ids` argument is ignored.
#' @param n_future_cv Number of time slices to hold out for future cross validation. Defaults
#' to floor(n_future_cv/3). Predictions are made for the last n_future_cv time steps
#' @param interval uncertainty interval for prediction. Can be one of "confidence" or "prediction",
#' and defaults to "confidence"
#' @param ... not used
#' @return A list object, containing `cv_pred` (a matrix of predictions), `cv_se` (a matrix of SEs),
#' `fold_ids` (a matrix of fold ids used as data), and `df` (a dataframe containing
#' the original data, predictions, SEs, and folds)
#'
#' @export
#' @examples
#' \donttest{
#' dat <- t(harborSealWA)
#' dat <- dat[2:4, ] # remove the year row
#' # fit a model with 1 hidden state and 3 observation time series
#' # cross validation here is random, 10 folds
#' fit <- MARSScv(dat, model = list(
#'   Z = matrix(1, 3, 1),
#'   R = "diagonal and equal"
#' ))
#'
#' # second, demonstrate passing in pre-specified folds
#' fold_ids <- matrix(
#'   sample(1:5, size = nrow(dat) * ncol(dat), replace = TRUE),
#'   nrow(dat), ncol(dat)
#' )
#' fit <- MARSScv(dat, model = list(
#'   Z = matrix(1, 3, 1),
#'   R = "diagonal and equal"
#' ), fold_ids = fold_ids)
#'
#' # third, illustrate future cross validation
#' fit <- MARSScv(dat, model = list(
#'   Z = matrix(1, 3, 1),
#'   R = "diagonal and equal"
#' ), future_cv = TRUE, n_future_cv = 5)
#' }
MARSScv <- function(y,
                    model = NULL,
                    inits = NULL,
                    method = "kem",
                    form = "marxss",
                    silent = FALSE,
                    control = NULL,
                    fun.kf = c("MARSSkfas", "MARSSkfss"),
                    fold_ids = NULL,
                    future_cv = FALSE,
                    n_future_cv = floor(ncol(y) / 3),
                    interval = "confidence",
                    ...) {
  if (is.null(fold_ids)) {
    # if fold_ids isn't entered, ids randomly
    n_folds <- 10
    fold_ids <- matrix(0, nrow(y), ncol(y))
    for (i in 1:nrow(fold_ids)) {
      for (j in 1:ncol(fold_ids)) fold_ids[i, j] <- sample(1:10, size = 1)
    }
  } else {
    n_folds <- max(fold_ids, na.rm = T)
    if (nrow(fold_ids) != nrow(y)) stop("Error: the number of rows of fold_ids must be the same as y")
    if (ncol(fold_ids) != ncol(y)) stop("Error: the number of rows of fold_ids must be the same as y")
  }

  pred <- matrix(NA, nrow(y), ncol(y))
  pred_se <- matrix(NA, nrow(y), ncol(y))
  row_indx <- matrix(rep(1:nrow(y), ncol(y)), nrow(y), ncol(y))
  col_indx <- matrix(sort(rep(1:ncol(y), nrow(y))), nrow(y), ncol(y))

  if (future_cv == TRUE) {
    # bookkeeping for the future cv
    n_folds <- n_future_cv
    fold_ids <- ncol(y) - col_indx + 1
  }

  for (k in 1:n_folds) {
    # drop data from this fold by making it NA
    if (future_cv == FALSE) {
      yfit <- y
      yfit[which(fold_ids == k)] <- NA
    } else {
      # drop out all data points
      yfit <- y
      yfit[, (ncol(y) - k + 1):ncol(y)] <- NA
    }
    # fit the model
    fit <- MARSS(
      y = yfit,
      model = model,
      inits = inits,
      method = method,
      form = form,
      silent = silent,
      control = control,
      fun.kf = fun.kf,
      ...
    )
    # save predictions and SEs
    pred_df <- predict(fit, interval = interval)$pred
    pred_mat <- matrix(pred_df$estimate, nrow = nrow(y), byrow = TRUE)
    pred_se_mat <- matrix(pred_df$se, nrow = nrow(y), byrow = TRUE)
    pred[which(fold_ids == k)] <- pred_mat[which(fold_ids == k)]
    pred_se[which(fold_ids == k)] <- pred_se_mat[which(fold_ids == k)]
  }

  # also return these as a df for plotting, etc.
  df <- pred_df[, c(".rownames", "t", "y")]
  df$cv_pred <- c(t(pred))
  df$cv_se <- c(t(pred_se))
  df$fold_id <- c(t(fold_ids))
  return(list(cv_pred = pred, cv_se = pred_se, fold_ids = fold_ids, df = df))
}
