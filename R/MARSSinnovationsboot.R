########################################################################################################################
#   MARSSinnovationsboot.r
#   This bootstraps multivariate time series data using Stoffer and Wall algorithm
#   It creates bootstrap data via sampling from the standardized innovations matrix
#   In the MARSS code, this is referred to as the nonparametric bootstrap.  Strictly speaking, it is not nonparametric.
########################################################################################################################
MARSSinnovationsboot = function(MLEobj, nboot=1000, minIndx=3 ) {
   if(any(is.na(MLEobj$marss$data))) 
      stop("Stopped in MARSSinnovationsboot() because this algorithm resamples from the innovations and doesn't allow missing values.\n", call.=FALSE)

   # The following need to be present in model: parameter estimates + Sigma, Kt, Innov; the latter 3 are part of $kf
   if(any(is.null(MLEobj[["par"]]), is.null(MLEobj[["marss"]][["data"]])))
      stop("Stopped in MARSSinnovationsboot(). The passed in marssMLE object is missing par or model$data elements.\n", call.=FALSE)
   if(is.null(MLEobj[["kf"]]) | !is.array(MLEobj[["kf"]][["Innov"]]) | !is.array(MLEobj[["kf"]][["Sigma"]])){ kf=MARSSkfss(MLEobj)
   }else{ kf=MLEobj[["kf"]] } 
   
   ## Rename things for code readability
   modelObj=MLEobj[["marss"]]
   par.dims=attr(modelObj,"model.dims")
   TT = par.dims[["data"]][2]
   m = par.dims[["x"]][1]   
   n = par.dims[["y"]][1]
   f=modelObj[["fixed"]]
   d=modelObj[["free"]]
   pari=parmat(MLEobj,t=1)

   #### make a list of time-varying parameters
   time.varying = list()
   for(elem in attr(modelObj,"par.names")) {
    if( (dim(d[[elem]])[3] == 1) & (dim(f[[elem]])[3] == 1)){
        time.varying[[elem]] = FALSE
      }else{ time.varying[[elem]] = TRUE }  #not time-varying
   }
      
   ##### Set holders for output
   newData = matrix(NA, n, TT)    #store as written for state-space eqn with time across columns
   newStates = matrix(NA, m, TT+1) #store as written for state-space eqn with time across columns  
   boot.data = array(NA, dim=c(par.dims[["data"]], nboot))  
   boot.states = array(NA, dim=c(par.dims[["x"]], nboot))  

    # Calculate the sqrt of sigma matrices, so they don't have to be computed 5000+ times
    sigma.Sqrt = array(0, dim=c(n, n, TT))
    BKS = array(0, dim=c(m, n, TT)) 
    for(i in 1:TT) {
      std.innovations = stdInnov(kf$Sigma, kf$Innov)  # standardized innovations; time across rows
      eig = eigen(kf$Sigma[,,i])   # calculate inv sqrt(sigma[1])
      sigma.Sqrt[,,i] = eig$vectors %*% makediag(sqrt(eig$values)) %*% t(eig$vectors)
      if(time.varying$B & i>1) pari$B = parmat(MLEobj,"B",t=i)$B
      BKS[,,i] = pari$B %*% kf$Kt[,,i] %*% sigma.Sqrt[,,i]   # this is m x n
    }

   for(i in 1:nboot){
      # Then the bootstrap algorithm proceeds
      # Stoffer & Wall suggest not sampling from innovations 1-3 (not much data)
      minIndx = ifelse(TT > 5, minIndx, 1)
      samp = sample(seq(minIndx+1, TT), size = (TT-minIndx), replace = TRUE)
      e = as.matrix(std.innovations)
      e[,1:minIndx] = as.matrix(std.innovations[,1:minIndx])
      e[,(minIndx+1):TT] = as.matrix(std.innovations[,samp])
      # newStates is a[] in the writeup by EH
      newStates[,1] = pari$x0
      #newStates[,2] = B %*% x0 + U + B %*% Kt[,,1] %*% sigma.Sqrt[,,1] %*% e[,1]
      for(t in 1:TT) {
         if(time.varying$B & i>1) pari$B = parmat(MLEobj,"B",t=i)$B
         if(time.varying$U & i>1) pari$U = parmat(MLEobj,"U",t=i)$U
         if(time.varying$Z & i>1) pari$Z = parmat(MLEobj,"Z",t=i)$Z
         if(time.varying$A & i>1) pari$A = parmat(MLEobj,"A",t=i)$A
         newStates[,t+1] = pari$B %*% newStates[,t,drop=FALSE] + pari$U + BKS[,,t] %*% e[,t,drop=FALSE]
         newData[,t] = pari$Z %*% newStates[,t,drop=FALSE] + pari$A + sigma.Sqrt[,,t] %*% e[,t,drop=FALSE]
      }

      newStates = newStates[,2:(TT+1)]
      boot.data[,,i] = newData
      boot.states[,,i] = newStates
      # reset newStates to its original dim
      newStates = matrix(NA, m, TT+1)
   }
     return(list(boot.states=boot.states, boot.data=boot.data, marss=modelObj, nboot=nboot))
}

######################################################################################################################
#   stdInnov
######################################################################################################################
stdInnov = function(SIGMA, INNOV) {
   # This function added by EW Nov 3, 2008
   # SIGMA is covariance matrix, E are original innovations
   TT = dim(INNOV)[2]
   n = dim(INNOV)[1]
   SI = matrix(0, n, TT)
   for(i in 1:TT) {
      a = SIGMA[,,i]
      a.eig = eigen(a)
      a.sqrt = a.eig$vectors %*% makediag((sqrt(a.eig$values))) %*% solve(a.eig$vectors)
      SI[,i] = chol2inv(chol(a.sqrt)) %*% INNOV[,i]  # eqn 1, p359 S&S
   }
   SI[which(INNOV == 0)]=0   # make sure things that should be 0 are 0
   return(SI)
}
