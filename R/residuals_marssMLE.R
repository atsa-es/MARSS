residuals.marssMLE = function(object,...){
  MLEobj=object
  #Reference page 9 in Messy Time Series
  #By definition there are no residuals for t=1
  model.dims=attr(MLEobj$marss,"model.dims")
  TT = model.dims[["x"]][2]
  m = model.dims[["x"]][1]
  n = model.dims[["y"]][1]
  rt = matrix(0,m,TT)
  ut = matrix(0,n,TT)
  et = st.et = matrix(0,n+m,TT)
  var.et = array(0,dim=c(n+m,n+m,TT))
  Nt = array(0,dim=c(m,m,TT))
  Mt = array(0,dim=c(n,n,TT))
  Jt = matrix(0,m,TT)
  #MARSSkfas doesn't output Innov, Sigma or Kt so might need to run MARSSkfss to get those
  if(is.null(MLEobj[["kf"]]) | is.null(MLEobj$kf$Innov) | is.null(MLEobj$kf$Sigma) | is.null(MLEobj$kf$Kt)) MLEobj$kf=MARSSkfss(MLEobj)
  #MARSSkfas sets these to a character warning, so not NULL; add this to catch that
  if(!is.array(MLEobj$kf$Innov) | is.array(MLEobj$kf$Sigma) | is.array(MLEobj$kf$Kt)) MLEobj$kf=MARSSkfss(MLEobj)
  vt = MLEobj$kf$Innov
  Ft = MLEobj$kf$Sigma
  for(t in seq(TT,2,-1)){
    H = matrix(0,m,n+m)
    Qt=parmat(MLEobj,"Q",t=t)$Q
    if(!any(takediag(Qt)==0)) H[,(n+1):(n+m)] = t(chol(Qt))
    G = matrix(0,n,n+m)
    Rt=parmat(MLEobj,"R",t=t)$R
    if(!any(takediag(Rt)==0)) G[,1:n] = t(chol(Rt))
    Tt = parmat(MLEobj,"B",t=t)$B
    Zt = parmat(MLEobj,"Z",t=t)$Z
    Ftinv = solve(Ft[,,t])
    Kt = Tt%*%matrix(MLEobj$kf$Kt[,,t],m,n) #R is dropping the dims so we force it to be nxm
    Lt = Tt - Kt%*%Zt
    Jt = H-Kt%*%G  
    ut[,t] = Ftinv%*%vt[,t,drop=FALSE]-t(Kt)%*%rt[,t,drop=FALSE]
    rt[,t-1] = t(Zt)%*%ut[,t,drop=FALSE]+t(Tt)%*%rt[,t,drop=FALSE]
    Mt[,,t] = Ftinv + t(Kt)%*%Nt[,,t]%*%Kt
    Nt[,,t-1] = t(Zt)%*%Ftinv%*%Zt + t(Lt)%*%Nt[,,t]%*%Lt
    et[,t] = t(G)%*%ut[,t,drop=FALSE] + t(H)%*%rt[,t,drop=FALSE]
    var.et[,,t] = t(G)%*%Ftinv%*%G + t(Jt)%*%Nt[,,t]%*%Jt
    st.er.et = matrix(0,n+m,n+m)
    if(!any(takediag(var.et[1:n,1:n,t])==0)) st.er.et[1:n,1:n] = solve(chol(var.et[1:n,1:n,t]))
    if(!any(takediag(var.et[(n+1):(n+m),(n+1):(n+m),t])==0)) st.er.et[(n+1):(n+m),(n+1):(n+m)] = solve(chol(var.et[(n+1):(n+m),(n+1):(n+m),t]))
    st.et[,t] = st.er.et%*%et[,t,drop=FALSE]
  }
  return(list(model.residuals=et[1:n,,drop=FALSE], state.residuals=et[(n+1):(n+m),,drop=FALSE], residuals=et, std.residuals=st.et, var.residuals=var.et))
}