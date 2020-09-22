accuracy.marssPredict <- function (f, x, test=NULL, type="ytt1", verbose=FALSE, ...) 
{
  out <- accuracy.marssMLE(f$model, test=test, type=type, verbose=verbose)
  
  if(!missing(x)){
    dotest <- TRUE
    if(f$h==0){
      message("Test data (x) passed in but object is not a forecast (h=0). Test data ignored.")
      dotest <- FALSE
    }
    if(dotest){ # f is a marssPredict object
      # set up the predictions
      n <- attr(f$model$model, "model.dims")[["y"]][1]
      h <- f$h

      TT <- attr(f$model$model, "model.dims")[["x"]][2]
      loc <- which(colnames(f$pred)=="se")
      testset <- subset(f$pred, f$t >TT)[,c(".rownames","y","estimate")]
      colnames(testset) <- c(".rownames","y",".fitted")
      testset$t <- rep(1:h,n)
      
      if(is.vector(x)) x <- matrix(x, nrow=1)
      if(inherits(x, "ts")) x <- t(x)
      if(!is.matrix(x) && !is.data.frame(x))
        stop("accuracy.marssMLE: Test data must be a matrix or data frame.\n", call.=FALSE)
      Y.names <- attr(f$model$model, "Y.names")
      if(is.matrix(x)){
        if(nrow(x) != n) stop(paste0("accuracy.marssPredict: Test data must have the same number of rows (n=",n,") as training data.\n"), call.=FALSE)
        if(ncol(x) != h) stop(paste0("accuracy.marssPredict: Test data must have the same columns as time steps (h=",h,") as in the forecast.\n"), call.=FALSE)
        if(is.null(rownames(x))){
          message("Training data rownames being used as rownames for test data.")
          rownames(x) <- attr(f$model$model, "Y.names")
        }
        val <- Y.names %in% rownames(x)
        if(!all(val)) stop(paste0("accuracy.marssPredict: Test data is missing ", paste(Y.names[!val], collapse=" ,"),".\n"), call.=FALSE)
        x <- x[match(Y.names, rownames(x)),,drop=FALSE] #match ordering if off
        x <- data.frame(
          .rownames=rep(rownames(x), each=h), 
          y=vec(t(x)), t=rep(1:h,n),
          stringsAsFactors = FALSE)
      }else{
        if(!all(c(".rownames", "y") %in% colnames(x))) stop("accuracy.marssPredict: If data frame, test data have rows .rownames and estimate.\n", call.=FALSE)
        val <- Y.names %in% x$.rownames
        if(!all(val)) stop(paste0("accuracy.marssPredict: Test data is missing ", paste(Y.names[!val], collapse=" ,"),".\n"), call.=FALSE)
        if(nrow(x) != h*n) stop(paste0("accuracy.marssPredict: Test data must have the same number of time steps (h=", h, ") as the forecast.\n"), call.=FALSE)
        if(!all(x$.rownames==testset$.rownames))
          loc <- match( paste(testset$.rownames, testset$t), paste(x$.rownames, x$t))
          x <- x[loc,]
      }
      testset$y <- x$y
      testout <- aout(testset, test)
      rout <- "Test set"
      rownames(testout) <- rout
      if(verbose){
        for(i in unique(testset$.rownames)){
          testout <- rbind(testout, aout(subset(testset, testset$.rownames==i), test))
          rout <- c(rout, paste(" ",i))
        }
        rownames(testout) <- rout
      }
      out <- rbind(out, testout)
    }
  }
  return(out)
}

accuracy.marssMLE <- function (f, x=NULL, test=NULL, type="ytt1", verbose=FALSE, ...)
{
  if(!missing(x)){
    message("Test data (x) passed in but object is not a forecast. x ignored.")
  }
  rout <- c("Training set")
  fx <- fitted(f, type=type)
  out <- aout(fx, test)
  rownames(out) <- rout
  n <- attr(f$model, "model.dims")[["y"]][1]
  if(!verbose && n!=1) return(out)
  for(i in unique(fx$.rownames)){
    out <- rbind(out, aout(subset(fx, fx$.rownames==i), test))
    rout <- c(rout, paste(" ",i))
  }
  rownames(out) <- rout
  return(out)
}

aout <- function(fx, test){ #fx is a data frame
  dx <- fx$y
  fits <- fx$.fitted
  res <- dx - fits
  if (is.null(test)) {
    test <- unique(fx$t)
  }
  if (!all(test %in% fx$t)) {
    stop("accuracy: test elements must be within sample.\n", call.=FALSE)
  }
  res <- res[fx$t %in% test]
  dx <- dx[fx$t %in% test]
  pe <- res/dx * 100
  me <- mean(res, na.rm = TRUE)
  mse <- mean(res^2, na.rm = TRUE)
  mae <- mean(abs(res), na.rm = TRUE)
  mape <- mean(abs(pe), na.rm = TRUE)
  mpe <- mean(pe, na.rm = TRUE)
  scale <- mean(abs(dx - mean(dx, na.rm = TRUE)), na.rm = TRUE)
  mase <- mean(abs(res/scale), na.rm = TRUE)
  
  if (length(res) > 1) {
    r1 <- acf(res, plot = FALSE, lag.max = 2, na.action = na.pass)$acf[2, 1, 1]
  }
  else {
    r1 <- NA
  }
  out <- matrix(c(me, sqrt(mse), mae, mpe, mape, mase, r1),nrow=1)
  colnames(out) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", "MASE", "ACF1")
  return(out)
}
