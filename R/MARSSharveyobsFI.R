####################################################################################
#   MARSSharveyobsFI function
#   Recursion to compute the observed Fisher Information matrix; Harvey (1989, pages 142-143)
#   With modification for missing values
#   Reference Holmes, E. E. (2014). Computation of standardized residuals for (MARSS) models. Technical Report. arXiv:1411.0045 [stat.ME]
####################################################################################
MARSSharveyobsFI <- function(MLEobj) {
  paramvector <- MARSSvectorizeparam(MLEobj)
  par.names <- names(paramvector)
  num.p <- length(paramvector)
  if (num.p == 0) {
    return(obsFI = matrix(0, 0, 0))
  } # no estimated parameters

  condition.limit <- 1E10
  condition.limit.Ft <- 1E5 # because the Ft is used to compute LL and LL drop limit is about 2E-8

  MODELobj <- MLEobj$marss
  n <- dim(MODELobj$data)[1]
  TT <- dim(MODELobj$data)[2]
  m <- dim(MODELobj$fixed$x0)[1]

  # create the YM matrix
  YM <- matrix(as.numeric(!is.na(MODELobj$data)), n, TT)
  # Make sure the missing vals in y are zeroed out if there are any
  y <- MODELobj$data
  y[YM == 0] <- 0

  if (MODELobj$tinitx == 1) {
    init.state <- "x10"
  } else {
    init.state <- "x00"
  }
  msg <- NULL
  # Construct needed identity matrices
  I.m <- diag(1, m)
  I.n <- diag(1, n)

  # Get the Kalman filter elements
  kf <- try(MARSSkfss(MLEobj), silent=TRUE) # kfss needed to get Sigma
  if( inherits(kf, "try-error")) stop("Stopped in MARSSharveyobsFI(). MARSSkfss does not run for this model. Try a numerical Hessian?", call.=FALSE)
  xtt <- kf$xtt
  xtt1 <- kf$xtt1
  Vtt <- kf$Vtt
  Vtt1 <- kf$Vtt1
  vt <- kf$Innov # these are innovations
  Ft <- kf$Sigma

  # initialize matrices
  # d values needed in recursion
  dxtt1 <- matrix(0, m, 1)
  dVtt1 <- matrix(0, m, m)
  dxt1t1 <- matrix(0, m, num.p)
  dVt1t1 <- array(0, dim = c(m, m, num.p))
  # need to store for each parameter since this is dvt/dtheta_i
  dvt <- matrix(0, n, num.p)
  dFt <- array(0, dim = c(n, n, num.p))
  obsFI <- matrix(0, num.p, num.p)
  rownames(obsFI) <- par.names
  colnames(obsFI) <- par.names

  ##############################################
  # Set up the parameters; Same code as in MARSSkf
  ##############################################
  # Note diff in param names B=T, u=c, a=d aka My notation (L)=Harvey notation (R)
  model.elem <- attr(MODELobj, "par.names")
  dims <- attr(MODELobj, "model.dims")
  time.varying <- c()
  for (el in model.elem) {
    if ((dim(MODELobj$free[[el]])[3] != 1) || (dim(MODELobj$fixed[[el]])[3] != 1)) {
      time.varying <- c(time.varying, el)
    }
  }
  pari <- parmat(MLEobj, t = 1)
  Z <- pari$Z
  A <- pari$A
  B <- pari$B
  U <- pari$U
  x0 <- pari$x0
  R <- tcrossprod(pari$H %*% pari$R, pari$H)
  Q <- tcrossprod(pari$G %*% pari$Q, pari$G)
  V0 <- tcrossprod(pari$L %*% pari$V0, pari$L)
  # set up the partial d par mats
  dpari <- dparmat(MLEobj, t = 1)
  # which par are time.varying
  time.varying.par <- time.varying[time.varying %in% names(dpari)]

  # set up an all zero dpar
  dpar0 <- list()
  for (el in model.elem) dpar0[[el]] <- matrix(0, dims[[el]][1], dims[[el]][2])

  ##############################################
  # RECURSION
  ##############################################

  for (t in 1:TT) {
    if (length(time.varying) != 0) {
      # update the time.varying ones
      pari[time.varying] <- parmat(MLEobj, time.varying, t = t)
      Z <- pari$Z
      A <- pari$A
      B <- pari$B
      U <- pari$U # update

      dpari[time.varying.par] <- dparmat(MLEobj, time.varying.par, t = t)
    }

    pcntr <- 0 # counter
    for (el in model.elem) {
      dp <- dpar0
      p <- length(MLEobj$par[[el]])
      if (p == 0) next # no estimated parameters for this parameter matrix
      for (ip in 1:p) {
        pcntr <- pcntr + 1
        # must ensure this is a matrix
        dp[[el]] <- array(dpari[[el]][, , ip], dim = dims[[el]][1:2]) # each is an array of c(dim(el),p)


        dHRH <- tcrossprod(dp$H %*% pari$R, pari$H) + tcrossprod(pari$H %*% dp$R, pari$H) + tcrossprod(pari$H %*% pari$R, dp$H)
        dGQG <- tcrossprod(dp$G %*% pari$Q, pari$G) + tcrossprod(pari$G %*% dp$Q, pari$G) + tcrossprod(pari$G %*% pari$Q, dp$G)
        dLV0L <- tcrossprod(dp$L %*% pari$V0, pari$L) + tcrossprod(pari$L %*% dp$V0, pari$L) + tcrossprod(pari$L %*% pari$V0, dp$L)

        # missing value modifications for Z, A, R per S&S2006 eq 6.78
        if (any(YM[, t] == 0)) {
          Mt <- I.n
          Mt[YM[, t] == 0, ] <- 0 # much faster than makediag(YM)
          I.2 <- I.n - Mt
          Zt <- Mt %*% Z # If Y missing, that row is 0 in Zt
          dp$Zt <- Mt %*% dp$Z
          dp$At <- Mt %*% dp$A
          # We need to deal with Ft when there are missing values
          # Normally Ft=1 on the diagonal for missing vals;
          # that way the LL calc works out
          # Ft=E(vt vt) is not defined when there are missing values
          # So we want Ft row i and col i to be 0 when the i-th obs is missing.
          # this is different than Ft, which we set to 1 on diag,
          # but the diag set to 1 doesn't matter since Vtt1 will have 0s in
          # in row i col i so will 0 that 1 out
          dHRHt <- Mt %*% dHRH %*% Mt
        } else {
          Zt <- Z
          dp$Zt <- dp$Z
          dHRHt <- dHRH
          dp$At <- dp$A
        }

        # t=1 treatment depends on how you define the initial condition.
        # Either as x at t=1 or x at t=0
        if (t == 1) {
          if (init.state == "x00") {
            # derivs of eqns 6.19 and 6.20 in S&S2006 with x0=x_0^0 and V0=V_0^0
            dxtt1 <- dp$B %*% x0 + B %*% dp$x0 + dp$U #
            dVtt1 <- tcrossprod(dp$B %*% V0, B) + tcrossprod(B %*% dp$V0, B) + tcrossprod(B %*% V0, dp$B) + dGQG
          }
          if (init.state == "x10") { # Ghahramani treatment of initial states
            # uses x10 and has no x00 (pi is defined as x10 at t=1);
            # See Holmes 2012.
            dxtt1 <- dp$x0
            dVtt1 <- dp$V0
          }
        } else { # t!=1
          # derivs of eqns 6.19 and 6.20 in S&S2006
          dxtt1 <- dp$B %*% xtt[, t - 1, drop = FALSE] + B %*% dxt1t1[, pcntr, drop = FALSE] + dp$U
          dVtt1 <- tcrossprod(dp$B %*% Vtt[, , t - 1], B) + tcrossprod(B %*% dVt1t1[, , pcntr], B) +
            tcrossprod(B %*% Vtt[, , t - 1], dp$B) + dGQG
        }
        if (m != 1) dVtt1 <- symm(dVtt1)
        # in general Vtt1 is not symmetric but here it is since Vtt and Q are

        # equations 3.4.71b and 3.4.73 in Harvey 1989; store for each p
        # the At, Zt etc, denotes that the missing vals mod vrs is used (does not denote time)
        dvt[, pcntr] <- -Zt %*% dxtt1 - dp$Zt %*% xtt1[, t, drop = FALSE] - dp$At
        dFt[, , pcntr] <- tcrossprod(dp$Zt %*% Vtt1[, , t], Zt) + tcrossprod(Zt %*% dVtt1, Zt) +
          tcrossprod(Zt %*% Vtt1[, , t], dp$Zt) + dHRHt # will always be matrix; Vtt1 is always square

        # get the inv of Ft
        if (n == 1) {
          Ftinv <- pcholinv(matrix(Ft[, , t], 1, 1))
        } else {
          Ftinv <- pcholinv(Ft[, , t]) # pcholinv deals with 0s on diagonal
          Ftinv <- symm(Ftinv) # to enforce symmetry after chol2inv call
        }

        # equation 3.4.74a in Harvey 1989; sets up dxtt[,t-1] needed for dxtt1
        dxt1t1[, pcntr] <- dxtt1 +
          tcrossprod(dVtt1, Zt) %*% Ftinv %*% vt[, t, drop = FALSE] +
          tcrossprod(Vtt1[, , t], dp$Zt) %*% Ftinv %*% vt[, t, drop = FALSE] -
          tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% dFt[, , pcntr] %*% Ftinv %*% vt[, t, drop = FALSE] +
          tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% dvt[, pcntr, drop = FALSE]

        # equation 3.4.74b in Harvey 1989; sets up dVtt[,t-1] needed for dVtt1
        dVt1t1[, , pcntr] <- dVtt1 -
          tcrossprod(dVtt1, Zt) %*% Ftinv %*% Zt %*% Vtt1[, , t] -
          tcrossprod(Vtt1[, , t], dp$Zt) %*% Ftinv %*% Zt %*% Vtt1[, , t] +
          tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% dFt[, , pcntr] %*% Ftinv %*% Zt %*% Vtt1[, , t] -
          tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% dp$Zt %*% Vtt1[, , t] -
          tcrossprod(Vtt1[, , t], Zt) %*% Ftinv %*% Zt %*% dVtt1
        if (m != 1) dVt1t1[, , pcntr] <- symm(dVt1t1[, , pcntr]) # to ensure its symetric
        # Check that cVt1t1 is zero-ed out like Vtt when R had 0s
      } # end of ip to p
    } # end of el in elem
    tmp <- matrix(0, num.p, num.p)
    for (i in 1:num.p) {
      for (j in i:num.p) {
        tmp[i, j] <- sum(diag(Ftinv %*% dFt[, , i] %*% Ftinv %*% dFt[, , j])) / 2 + t(dvt[, i, drop = FALSE]) %*% Ftinv %*% dvt[, j, drop = FALSE]
        tmp[j, i] <- tmp[i, j]
      }
    }
    obsFI <- obsFI + tmp
  } # End of the recursion (for i to 1:TT)

  return(obsFI)
}

symm <- function(x) {
  t.x <- matrix(x, dim(x)[2], dim(x)[1], byrow = TRUE)
  x <- (x + t.x) / 2
  x
}

dparmat <- function(MLEobj, elem = c("B", "U", "Q", "Z", "A", "R", "x0", "V0", "G", "H", "L"), t = 1) {
  # returns a list where each el in elem is an element.  Returns a 3D matrix; one for each estimated parameter in the matrix
  # Only returns values for parameters with estimated elements
  # needs MLEobj$marss and MLEobj$par
  # f=MLEobj$marss$fixed
  # delem is the elem that theta_i is from.  i is from 1 to length(MLEobj$par$delem)
  model.loc <- "marss"
  model <- MLEobj[[model.loc]]
  pars <- MLEobj[["par"]]
  d <- model[["free"]]
  num.p <- length(MARSSvectorizeparam(MLEobj))
  par.mat <- list()
  dims <- attr(model, "model.dims")
  # Will be NULL unless the el has parameter elements being estimated
  if (num.p != 0) { # Something is estimated
    for (el in elem) { # this is for the list
      par.el <- pars[[el]] # estimated parameters for this el
      p <- length(par.el)
      if (p == 0) next # next el since this one has no estimated parameters
      par.mat[[el]] <- array(0, dim = c(dims[[el]][1:2], p))
      for (ip in 1:p) {
        dpar <- matrix(0, p, 1)
        dpar[ip] <- 1
        if (dim(d[[el]])[3] == 1) { # non-time-varying
          delem <- d[[el]]
        } else {
          delem <- d[[el]][, , t]
        }
        attr(delem, "dim") <- attr(delem, "dim")[1:2]
        par.mat[[el]][, , ip] <- matrix(delem %*% dpar, dims[[el]][1], dims[[el]][2])
      }
    }
  }
  return(par.mat)
}
