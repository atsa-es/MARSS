#######################################################################################################
#   MARSShatyt function
#   Expectations involving hatyt
#######################################################################################################
MARSShatyt = function( MLEobj ) {
  MODELobj = MLEobj[["marss"]]
  if(!is.null(MLEobj[["kf"]])){ kfList = MLEobj$kf
  }else{ kfList=MARSSkf(MLEobj) }
  model.dims=attr(MODELobj,"model.dims")
  n=model.dims$data[1]; TT=model.dims$data[2]; m=model.dims$x[1]

  y=MODELobj$data
  #create the YM matrix; need to keep track of location of missing
  #YM=matrix(as.numeric(!is.na(y)),n,TT) # TRUE is present ; FALSE is missing
  YM=!is.na(y)
  #Make sure the missing vals in y are zeroed out if there are any
  y[!YM]=0
  
  #set-up matrices for hatxt for 1:TT and 1:TT-1
  IIz=list(); IIz$V0=makediag(as.numeric(takediag(parmat(MLEobj, "V0", t=1)$V0)==0),m)
  #bad notation; should be hatxtT
  hatxt=kfList$xtT      #1:TT
  hatxtt1=kfList$xtt1      #1:TT
  E.x0 = (diag(1,m)-IIz$V0)%*%kfList$x0T+IIz$V0%*%parmat(MLEobj, "x0", t=1)$x0  #0:T-1
  hatxt1=cbind(E.x0,kfList$xtT[,1:(TT-1),drop=FALSE])
  hatxtp=cbind(kfList$xtT[,2:TT,drop=FALSE],NA)
  hatVt=kfList$VtT
  hatVtt1=kfList$Vtt1T
  #hatVtpt = array(NA,dim=dim(kfList$Vtt1T))
  #hatVtpt[,,1:(TT-1)] = kfList$Vtt1T[,,2:TT,drop=FALSE]
  
  msg=NULL
  
  #Construct needed identity matrices
  I.n = diag(1,n)
  
  #Note diff in param names from S&S;B=Phi, Z=A, A not in S&S
  time.varying = c(); pari=list()
  for(elem in c("R","Z","A")){  #only params needed for this function
    if( model.dims[[elem]][3] == 1 ){  #not time-varying
      pari[[elem]]=parmat(MLEobj, elem, t=1)[[elem]]
      if(elem=="R"){ 
        if(length(pari$R)==1) diag.R=unname(pari$R) else diag.R = takediag(unname(pari$R))
        #isDiagonal is rather expensive; this test is faster
        is.R.diagonal = all(pari$R[!diag(nrow(pari$R))] == 0) # = isDiagonal(pari$R)
      }
    }else{ time.varying = c(time.varying, elem) } #which elements are time varying
  }#end for loop over elem
  
  
  #initialize - these are for the forward, Kalman, filter
  # for notation purposes, 't' represents current point in time, 'TT' represents the length of the series
  #notation bad.  should be hatytT
  hatyt = hatytt1 = y # faster than matrix(0,n,TT)     
  hatOt = array(0,dim=c(n,n,TT))     
  hatyxt = hatyxtt1 = hatyxttp = array(0,dim=c(n,m,TT))
  
  for (t in 1:TT) {
    for(elem in time.varying){
      pari[[elem]]=parmat(MLEobj, elem, t=t)[[elem]] 
      if(elem=="R"){ 
        if(length(pari$R)==1) diag.R=unname(pari$R) else diag.R = takediag(unname(pari$R))
        #isDiagonal is rather expensive; this test is faster
        is.R.diagonal = all(pari$R[!diag(nrow(pari$R))] == 0)
        #is.R.diagonal = isDiagonal(pari$R)
      }
    }
    YMt=YM[,t]
    if(all(YMt)){  #none missing
      hatyt[,t]=y[,t,drop=FALSE]
      hatytt1[,t]=pari$Z%*%hatxtt1[,t,drop=FALSE]+pari$A
      hatOt[,,t]=base::tcrossprod(hatyt[,t,drop=FALSE]) 
#      hatOt[,,t]=hatyt[,t,drop=FALSE]%*%t(hatyt[,t,drop=FALSE],1,n)
      hatyxt[,,t]=base::tcrossprod(hatyt[,t,drop=FALSE], hatxt[,t,drop=FALSE])
#      hatyxt[,,t]=hatyt[,t,drop=FALSE]%*%t(hatxt[,t,drop=FALSE],1,m)
      hatyxtt1[,,t]=base::tcrossprod(hatyt[,t,drop=FALSE], hatxt1[,t,drop=FALSE])
#      hatyxtt1[,,t]=hatyt[,t,drop=FALSE]%*%t(hatxt1[,t,drop=FALSE],1,m)
      hatyxttp[,,t]=base::tcrossprod(hatyt[,t,drop=FALSE], hatxtp[,t,drop=FALSE])
    }else{
      I.2 = I.r = I.n; 
      I.2[YMt,]=0 #1 if YM=0 and 0 if YM=1
      I.r[!YMt | diag.R==0,]=0   #if Y missing or R = 0, then 0
      Delta.r=I.n
      if(is.R.diagonal) Delta.r = I.n-I.r
      if(!is.R.diagonal && any(YMt & diag.R!=0)){           
        mho.r = I.r[YMt & diag.R!=0,,drop=FALSE]
        t.mho.r = I.r[,YMt & diag.R!=0,drop=FALSE]
        Rinv = try(chol(mho.r%*%pari$R%*%t.mho.r))
        #Catch errors before entering chol2inv
        if(class(Rinv)=="try-error") {
          return(list(ok=FALSE, errors="Stopped in MARSShatyt: chol(R) error.\n" ) )      
        }
        Rinv=chol2inv(Rinv) 
        Delta.r = I.n- pari$R%*%t.mho.r%*%Rinv%*%mho.r
      }
      hatyt[,t]=y[,t,drop=FALSE] - Delta.r%*%(y[,t,drop=FALSE]-pari$Z%*%hatxt[,t,drop=FALSE]-pari$A)
      hatytt1[,t]=pari$Z%*%hatxtt1[,t,drop=FALSE]+pari$A
      #t.DZ = matrix(Delta.r%*%pari$Z,m,n,byrow=TRUE)
      DZ = Delta.r%*%pari$Z
      hatOt[,,t]=I.2%*%(Delta.r%*%pari$R + DZ%*%base::tcrossprod(hatVt[,,t], DZ))%*%I.2 + tcrossprod(hatyt[,t,drop=FALSE])
#      hatOt[,,t]=I.2%*%(Delta.r%*%pari$R+Delta.r%*%pari$Z%*%hatVt[,,t]%*%t.DZ)%*%I.2 + hatyt[,t,drop=FALSE]%*%t(hatyt[,t,drop=FALSE],1,n)
      hatyxt[,,t]=base::tcrossprod(hatyt[,t,drop=FALSE], hatxt[,t,drop=FALSE])+DZ%*%hatVt[,,t]
#      hatyxt[,,t]=hatyt[,t,drop=FALSE]%*%t(hatxt[,t,drop=FALSE],1,m)+Delta.r%*%pari$Z%*%hatVt[,,t]
      hatyxtt1[,,t]=base::tcrossprod(hatyt[,t,drop=FALSE],hatxt1[,t,drop=FALSE])+DZ%*%hatVtt1[,,t]
#      hatyxtt1[,,t]=hatyt[,t,drop=FALSE]%*%t(hatxt1[,t,drop=FALSE],1,m)+Delta.r%*%pari$Z%*%hatVtt1[,,t]
      if(t==TT){ hatyxttp[,,t] = NA }else{
      hatyxttp[,,t]=base::tcrossprod(hatyt[,t,drop=FALSE],hatxtp[,t,drop=FALSE]) + base::tcrossprod(DZ, Vtt1T[,,t+1])
      }
      #hatyxttp[,,t]=base::tcrossprod(hatyt[,t,drop=FALSE],hatxtp[,t,drop=FALSE])+Delta.r%*%pari$Z%*%t(hatVtpt[,,t])
    }
  } #for loop over time
  rtn.list=list(ytT = hatyt, OtT = hatOt, yxtT=hatyxt, yxt1T=hatyxtt1, yxttpT = hatyxttp, ytt1 = hatytt1) 
  return(c(rtn.list,list(ok=TRUE, errors = msg)))
}
