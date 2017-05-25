#######################################################################################################
#   MARSSkf function
#   Utility function to choose the Kalman filter and smoother
#######################################################################################################
MARSSkf = function( MLEobj, only.logLik=FALSE, return.lag.one=TRUE, return.kfas.model=FALSE ) {
if(is.null(MLEobj$par))
  stop("Stopped in MARSSkf(): par element of marssMLE object is required.\n")
#full=TRUE means to return all the kf info, this is currently done by MARSSkfss
if(MLEobj$fun.kf=="MARSSkfss") 
  return(MARSSkfss(MLEobj))
if(MLEobj$fun.kf=="MARSSkfas") 
  return( MARSSkfas(MLEobj, only.logLik=only.logLik, return.lag.one=return.lag.one, return.kfas.model=return.kfas.model ) )
return(list(ok=FALSE, errors="kf.function does not specify a valid Kalman filter and smoother function."))
}