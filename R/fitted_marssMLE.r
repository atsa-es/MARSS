###############################################################################################################################################
#  fitted method for class marssMLE. 
#  returns the fitted value of y conditioned on all the data or data up to t-1 if one.step.ahead=TRUE
##############################################################################################################################################
fitted.marssMLE <- function (object, ..., one.step.ahead=FALSE) {
  MLEobj=object
  if(is.null(MLEobj$par))
    stop("The marssMLE object does not have the par element.  Most likely the model has not been fit.")
  
  model.dims=attr(MLEobj$marss,"model.dims")
  TT = model.dims[["x"]][2]
  n = model.dims[["y"]][1]
  
  if(!one.step.ahead){ 
    hatxt = MLEobj$states
  }else{
    hatxt = MARSSkf(MLEobj)$xtt1
  }
  
  val = matrix(NA, n, TT)
  rownames(val) = attr(kemz.2$marss,"Y.names")
  
  for(t in 1:TT){
    Zt=parmat(MLEobj,"Z",t=t)$Z
    At=parmat(MLEobj,"A",t=t)$A    
    val[,t] = Zt%*%hatxt[,t,drop=FALSE]+At
  }
  
val
  
}  #end of fitted.marssMLE