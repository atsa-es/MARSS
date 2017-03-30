#######################################################################################################
#   MARSShessian functions
#   Adds Hessian, parameter var-cov matrix, and parameter mean to a marssMLE object
#   Returns par and fixed and free in chol form
#   The func MARSShessian.backtrans (below) transforms that to the original form
#   Note set up to allow unconstrained var-cov matrices but I don't think it really works for that
#######################################################################################################
MARSShessian.bad = function(MLEobj) {
  
  paramvector = MARSSvectorizeparam(MLEobj)
  
  kfNLL=function(paramvec, MLEobj){
    tmp.MLEobj = MARSSvectorizeparam(MLEobj, parvec=paramvec)
    negLL = MARSSkf(tmp.MLEobj, only.logLik=TRUE, return.lag.one=FALSE )$logLik
    
    -1*negLL
    #return(-1*logLik(tmp.MLEobj))
  }
  
  #Hessian and gradient
  emhess = fdHess(paramvector, function(paramvector, MLEobj) kfNLL(paramvector, MLEobj), MLEobj)
  MLEobj$Hessian = emhess$Hessian
  rownames(MLEobj$Hessian)=names(paramvector)
  colnames(MLEobj$Hessian)=names(paramvector)
  MLEobj$gradient = emhess$gradient
  
  parSigma = try(solve(MLEobj$Hessian), silent=TRUE)
  if(inherits(parSigma, "try-error")) {
    warning("MARSShessian: Hessian could not be inverted to compute the parameter var-cov matrix")
    parSigma=NULL
  }
  MLEobj$parSigma = parSigma
  MLEobj$parMean = paramvector
  
  #This is the UNTRANSFORMED MLEobj
  return(MLEobj)
}

