residuals.marssMLE = function(object,..., Harvey=FALSE, normalize=FALSE){
  #These are the residuals and their variance conditioned on all the data
  #Harvey=TRUE uses Harvey et al (1998) algorithm to compute these
  #Harvey=FALSE uses the straight smoother output
  # model.residuals y(t|yT)-Zx(t|yT)-a
  # state.residuals x(t|yT)-Bx(t-1|yT)-u
  # residuals 
  # var.residuals variance of above over Y
  #    except for missing values, var for y_i missing should be 0 but if Harvey=TRUE, returns R(i,i)
  # std.residuals psolve(chol(var.residuals))%*%model.residuals
  # normalize means to divide residuals by H%*%t(chol(R)) or G%*%t(chol(Q)) so they are scaled to I
  # Note, I think there is a problem with the Harvey algorithm when the variance of the state residuals (Q)
  # is non-diagonal and there are missing values; it can become non-invertable
  MLEobj=object
  model.dims=attr(MLEobj$marss,"model.dims")
  TT = model.dims[["x"]][2]
  m = model.dims[["x"]][1]
  n = model.dims[["y"]][1]
  y=MLEobj$marss$data
  #set up holders
  et = st.et = matrix(0,n+m,TT)
  var.et = array(0,dim=c(n+m,n+m,TT))
  
  #### make a list of time-varying parameters
  time.varying = list()
  for(elem in attr(MLEobj[["marss"]],"par.names")) {
    if( model.dims[[elem]][3]==1 ){
      time.varying[[elem]] = FALSE
    }else{ time.varying[[elem]] = TRUE }
  }
  
  #MARSSkfas doesn't output Innov, Sigma or Kt so might need to run MARSSkfss to get those
  if(is.null(MLEobj[["kf"]]) || is.null(MLEobj$kf$Innov) || is.null(MLEobj$kf$Sigma) || is.null(MLEobj$kf$Kt)) 
    kf=MARSSkfss(MLEobj)
  #MARSSkfas sets these to a character warning, so not NULL; add this to catch that
  if(!is.array(MLEobj$kf$Innov) || !is.array(MLEobj$kf$Sigma) || !is.array(MLEobj$kf$Kt)) 
    kf=MARSSkfss(MLEobj)
  Ey = MARSShatyt(MLEobj)
  
  #For debugging purposes I leave in the Harvey et al 1998 algorithm but default is to not use that
  if(!Harvey){
    #model.et will be 0 where no data E(y)-modeled(y)
    model.et = Ey$ytT-fitted(MLEobj, type="y", conditioning="T") #model residuals
    et[1:n,]=model.et
    
    for(t in 1:TT){
      #model residuals
      Rt=parmat(MLEobj,"R",t=t)$R #returns matrix
      Ht=parmat(MLEobj,"H",t=t)$H
      #varcov with H and R togethter
      Rt=Ht%*%Rt%*%t(Ht)      
      
      Zt=parmat(MLEobj,"Z",t=t)$Z
      
      
      #model.et defined outside for loop
      
      #compute the variance of the residuals and state.et
      St = Ey$yxtT[,,t] - tcrossprod(Ey$ytT[,t,drop=FALSE],kf$xtT[,t,drop=FALSE])
      tmpvar.et = Rt-Zt%*%kf$VtT[,,t]%*%t(Zt) + St%*%t(Zt) + Zt%*%t(St)
      
      if(t<TT){ #fill in var.et for t (model resid)
        Qtp=parmat(MLEobj,"Q",t=t+1)$Q
        Gtp=parmat(MLEobj,"G",t=t+1)$G
        Qtp=Gtp%*%Qtp%*%t(Gtp)
        
        Btp=parmat(MLEobj,"B",t=t+1)$B
        utp=parmat(MLEobj,"U",t=t+1)$U

        Sttp = Ey$yxttpT[,,t] - tcrossprod(Ey$ytT[,t,drop=FALSE],kf$xtT[,t+1,drop=FALSE])
        cov.et = Zt%*%t(kf$Vtt1T[,,t+1]) - Zt%*%kf$VtT[,,t]%*%t(Btp) - Sttp + St%*%t(Btp)
        tmpvar.state.et = Qtp-kf$VtT[,,t+1]-Btp%*%kf$VtT[,,t]%*%t(Btp)+kf$Vtt1T[,,t+1]%*%t(Btp)+Btp%*%t(kf$Vtt1T[,,t+1])
        
        et[(n+1):(n+m),t]=kf$xtT[,t+1]-Btp%*%kf$xtT[,t]-utp
      }else{
        cov.et = matrix(0,n,m)
        tmpvar.state.et = matrix(0,m,m)
      }
      var.et[1:n,,t] = cbind(tmpvar.et, cov.et)
      var.et[(n+1):(n+m),,t] = cbind(t(cov.et),tmpvar.state.et)
      
      if(normalize){
        Qpinv = matrix(0,m,n+m)
        Rinv = matrix(0,n,n+m)
        if(t<TT){
          Qtp=parmat(MLEobj,"Q",t=t+1)$Q
          Qpinv[,(n+1):(n+m)]=psolve(t(pchol(Qtp)))
        }
        Rinv[,1:n] = psolve(t(pchol(Rt)))
        RQinv=rbind(Rinv,Qpinv) #block diag matrix
        et[,t]=RQinv%*%et[,t]
        var.et[,,t] = RQinv%*%var.et[,,t]%*%t(RQinv)
      }
      
    }
    
  }else{ #use Harvey algorithm
    #Reference page 112-133 in Messy Time Series
    #Reference de Jong and Penzer 1998; with model transformed so sigma^2 = 1
    #refs in man file
    #NOTATION here uses that in Messy Time Series (Harvey et al 1998) with some differences
    #MARSS uses Koopman's terminology where w_t = G%*%w, Harvey uses H%*%w
    #By definition residual = 0 at time t, where x_0 is defined at t=0 or t=1, **when x_0 is estimated**
    #because x_0 is estimated and the max L will be when residual is 0
    #x_TT = Bx_{TT-1}+u_{TT-1}+w_{TT-1}, so w_TT cannot be computed (you'd need X_{TT+1}
    #so state residual at t=TT is NA
    #The state residuals are diff(states) but scaled by t(chol(Q))
    rt = matrix(0,m,TT)
    ut = matrix(0,n,TT)
    Nt = array(0,dim=c(m,m,TT))
    Mt = array(0,dim=c(n,n,TT))
    Jt = matrix(0,m,TT)
    vt = kf$Innov
    # If they are time-varying, Q, G and B at t=T will not appear (cancelled out by r_T and N_T = 0). 
    # Set for t=1. Will update in for loop if time-varying.
    pari = parmat(MLEobj,t=1)
    Qtp = pari[["Q"]]
    Gtp = pari[["G"]]
    Ttp = pari[["B"]]
    Zt = pari[["Z"]]
    Rt = pari[["R"]]
    Ht = pari[["H"]]
    for(t in seq(TT,1,-1)){
      #define all the, potential time-varying, parameters
      if(t<TT){
        if(time.varying[["Q"]]) Qtp=parmat(MLEobj,"Q",t=t+1)$Q
        if(time.varying[["G"]]) Gtp=parmat(MLEobj,"G",t=t+1)$G
        if(time.varying[["B"]]) Ttp = parmat(MLEobj,"B",t=t+1)$B
      }
      if(time.varying[["Z"]]) Zt = parmat(MLEobj,"Z",t=t)$Z
      if(time.varying[["R"]]) Rt=parmat(MLEobj,"R",t=t)$R
      if(time.varying[["H"]]) Ht=parmat(MLEobj,"H",t=t)$H

      #implement missing values modifications per Shumway and Stoffer
      diag.Rt=diag(Rt)
      Rt[is.na(y[,t]),]=0
      Rt[,is.na(y[,t])]=0
      #diag(Rt)=diag.Rt
      Zt[is.na(y[,t]),]=0
      
      #create the m x n+m and n x n+m matrices
      Rstar = matrix(0,n,n+m)      
      Rstar[,1:n] = Ht%*%Rt%*%t(Ht)
      if(normalize){ 
        Rstar[,1:n] = Ht%*%t(pchol(Rt))
      }else{
        Rstar[,1:n] = Ht%*%Rt%*%t(Ht)
      }
      Qpstar = matrix(0,m,n+m)
      #MARSS uses Koopman's terminology where w_t = G%*%w, Harvey uses H%*%w
      if(normalize){ 
        Qpstar[,(n+1):(n+m)] = Gtp%*%t(pchol(Qtp))
      }else{
        Qpstar[,(n+1):(n+m)] = Gtp%*%Qtp%*%t(Gtp)
      }
      
      #Don't use kf$Sigma since that has (i,i)=1
      Ftinv = psolve(Zt%*%kf$Vtt1[,,t]%*%t(Zt)+Rt)
      
      #Harvey algorithm modified to return non-normalized errors
      Kt = Ttp%*%matrix(kf$Kt[,,t],m,n) #R is dropping the dims so we force it to be mxn
      Lt = Ttp - Kt%*%Zt
      Jt = Qpstar-Kt%*%Rstar  
      ut[,t] = Ftinv%*%vt[,t,drop=FALSE]-t(Kt)%*%rt[,t,drop=FALSE]
      if(t>1) rt[,t-1] = t(Zt)%*%ut[,t,drop=FALSE]+t(Ttp)%*%rt[,t,drop=FALSE]
      #Mt[,,t] = Ftinv + t(Kt)%*%Nt[,,t]%*%Kt #not used
      if(t>1) Nt[,,t-1] = t(Zt)%*%Ftinv%*%Zt + t(Lt)%*%Nt[,,t]%*%Lt
      et[,t] = t(Rstar)%*%ut[,t,drop=FALSE] + t(Qpstar)%*%rt[,t,drop=FALSE]
      #see deJong and Penzer 1998, page 800, right column, halfway down
      var.et[,,t] = t(Rstar)%*%Ftinv%*%Rstar + t(Jt)%*%Nt[,,t]%*%Jt
      
    }
  }
  
  #prepare standardized residuals
  for(t in 1:(TT-1)){
    tmpvar = sub3D(var.et, t=t)
    resids = et[,t,drop=FALSE]
    #don't includ values for resids if there is no residual (no data)
    is.miss = c(is.na(y[,t]),rep(FALSE,m))
    resids[is.miss]=0

    tmpvar[abs(tmpvar)<sqrt(.Machine$double.eps)]=0
    
    #psolve and pchol deal with 0s on diagonal
    #wrapped in try to prevent crashing if inversion not possible
    tmpchol=try(pchol(tmpvar), silent=TRUE)
    if( inherits(tmpchol, "try-error") ){ 
      st.et[,t]=NA
      cat(paste("warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =",t,". See MARSSinfo(\"residvarinv\")\n"))
      next
    }
    tmpcholinv=try(psolve(tmpchol), silent=TRUE)
    if(inherits(tmpcholinv, "try-error")){
      st.et[,t]=NA
      cat(paste("warning: the variance of the residuals at t =", t, "is not invertible.  NAs returned for std.residuals at t =",t,"\n"))
      next
    }
    st.et[,t] = tmpcholinv%*%resids
    st.et[is.miss,t]=NA
  }
  
  # the state.residual at the last time step is NA because it is x(T+1) - f(x(T)) and T+1 does not exist.  For the same reason, the var.residuals at TT will have NAs
  et[(n+1):(n+m),TT] = NA
  var.et[,(n+1):(n+m),TT] = NA
  var.et[(n+1):(n+m),,TT] = NA
  st.et[,TT] = NA
  
  # add rownames
  Y.names = attr(MLEobj$model,"Y.names")
  X.names = attr(MLEobj$model,"X.names")
  rownames(et)=rownames(st.et)=rownames(var.et)=colnames(var.et)=c(Y.names, X.names)
  
  return(list(model.residuals=et[1:n,,drop=FALSE], state.residuals=et[(n+1):(n+m),,drop=FALSE], residuals=et, std.residuals=st.et, var.residuals=var.et))
}