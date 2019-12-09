residuals.t1.marssMLE = function(object,..., Harvey=FALSE, normalize=FALSE){
# residuals of one.step.ahead preditions
  # model.residuals y(t|yt1)-Zx(t|yt1)-a
  # state.residuals x(t|yt1)-Bx(t-1|yt1)-u
  # residuals 
  # var.residuals variance of above over Y
  #    except for missing values, var for y_i missing should be 0 but if Harvey=TRUE, returns R(i,i)
  # std.residuals psolve(chol(var.residuals))%*%model.residuals
  # normalize means to divide residuals by H%*%t(chol(R)) or G%*%t(chol(Q)) so they are scaled to I

  MLEobj=object
  model.dims=attr(MLEobj$marss,"model.dims")
  TT = model.dims[["x"]][2]
  m = model.dims[["x"]][1]
  n = model.dims[["y"]][1]
  y=MLEobj$marss$data
  #set up holders
  et = st.et = matrix(0,n+m,TT)
  var.et = array(0,dim=c(n+m,n+m,TT))
  
  #MARSSkfas doesn't output Innov, Sigma or Kt so might need to run MARSSkfss to get those
  if(is.null(MLEobj[["kf"]]) || is.null(MLEobj$kf$Innov) || is.null(MLEobj$kf$Sigma) || is.null(MLEobj$kf$Kt)) 
    kf=MARSSkfss(MLEobj)
  #MARSSkfas sets these to a character warning, so not NULL; add this to catch that
  if(!is.array(MLEobj$kf$Innov) || !is.array(MLEobj$kf$Sigma) || !is.array(MLEobj$kf$Kt)) 
    kf=MARSSkfss(MLEobj)
  Ey = MARSShatyt(MLEobj)
  
    #model.et will be 0 where no data since E(y)=modeled(y)
    #otherwise it is y-modeled(y)
    model.et = y-fitted(MLEobj, conditioning="t-1") #model residuals
    model.et[is.na(y)]=0
    et[1:n,]=model.et
    
    for(t in 1:TT){
      #model residuals
      Rt=parmat(MLEobj,"R",t=t)$R #returns matrix
      Ht=parmat(MLEobj,"H",t=t)$H
      #varcov with H and R togethter
      Rt=Ht%*%Rt%*%t(Ht)      
      Zt=parmat(MLEobj,"Z",t=t)$Z
      
      #model.et defined outside for loop
      
      #compute the variance of the model residuals
       tmpvar.model.et = Rt+Zt%*%kf$Vtt1[,,t]%*%t(Zt)
      
      if(t<TT){ #fill in var.et for state residual
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
      var.et[1:n,,t] = cbind(tmpvar.model.et, cov.et)
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