###############################################################################################################################################
#  fitted method for class marssMLE. 
#  returns the fitted value of y conditioned on all the data or data up to t-1 if one.step.ahead=TRUE
##############################################################################################################################################
fitted.marssMLE <- function (object, ..., one.step.ahead=FALSE, type=c("observations", "states")) {
  MLEobj=object
  if(is.null(MLEobj$par))
    stop("The marssMLE object does not have the par element.  Most likely the model has not been fit.")
  
  if(type=="observations"){
  model.dims=attr(MLEobj$marss,"model.dims")
  TT = model.dims[["x"]][2]
  n = model.dims[["y"]][1]
  
  if(!one.step.ahead){ 
    hatxt = MLEobj$states
  }else{
    hatxt = MARSSkf(MLEobj)$xtt1
  }
  
  val = matrix(NA, n, TT)
  rownames(val) = attr(MLEobj$marss,"Y.names")
  
  for(t in 1:TT){
    Zt=parmat(MLEobj,"Z",t=t)$Z
    At=parmat(MLEobj,"A",t=t)$A    
    val[,t] = Zt%*%hatxt[,t,drop=FALSE]+At
  }
  }
  if(type=="states"){
    model.dims=attr(MLEobj$marss,"model.dims")
    TT = model.dims[["x"]][2]
    n = model.dims[["x"]][1]
    
    if(!one.step.ahead){ 
      hatxt = MLEobj$states
    }else{
      hatxt = MARSSkf(MLEobj)$xtt1
    }
    
    val = matrix(NA, n, TT)
    rownames(val) = attr(MLEobj$marss,"X.names")
    
    x0 = coef(MLEobj, type="matrix")$x0
    Bt=parmat(MLEobj,"B",t=1)$B
    Ut=parmat(MLEobj,"U",t=1)$U    
    val[,1] = Bt%*%x0+Ut
    for(t in 2:TT){
      Bt=parmat(MLEobj,"B",t=t)$B
      Ut=parmat(MLEobj,"U",t=t)$U    
      val[,t] = Bt%*%hatxt[,t-1,drop=FALSE]+Ut
    }
  }
  
val
  
}  #end of fitted.marssMLE