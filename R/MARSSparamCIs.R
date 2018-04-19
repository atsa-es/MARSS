####################################################################################
#   MARSSparamCIs function
#   This returns CIs for ML parameter estimates
#   If method='hessian', uses either Harvey1989 (analytical) or fdHess (numerical) or optim (numerical)
####################################################################################
MARSSparamCIs = function(MLEobj, method="hessian", alpha=0.05, nboot=1000, silent=TRUE, hessian.fun="Harvey1989") {
  #this function expects a marssMLE object
  #it will add standard errors, biases, low/up CIs to MLEobj
  if(class(MLEobj)[1]!="marssMLE")
    stop("Stopped in MARSSparamCIs(). This function needs a marssMLE object.\n", call.=FALSE)

    if(!(method %in% c("hessian","parametric","innovations"))) stop("Stopped in MARSSparamCIs(). Current methods are hessian, parametric, innovations.\n", call.=FALSE)
  if(!(hessian.fun %in% c("Harvey1989","fdHess","optim"))) stop("Stopped in MARSSparamCIs(). Available functions for computing the Hessian are Harvey1989, fdHess, and optim.\n", call.=FALSE)  
  
  
  paramvec = MARSSvectorizeparam(MLEobj)
  if(length(paramvec)==0) stop("Stopped in MARSSparamCIs(). No estimated parameter elements.\n", call.=FALSE)
  
  paramnames=names(paramvec)
  
  #default
  if(method=="hessian"){
    MLEobj = MARSShessian(MLEobj, method=hessian.fun)
    
    vector.par.se =rep(NA,length(paramvec))
    if(is.null(MLEobj$parSigma)) {
      warning("MARSSparamCIs: No parSigma element returned by Hessian function.  See marssMLE object errors (MLEobj$errors)")
    }else{
      is.neg=diag(MLEobj$parSigma)<0
      if(any(is.neg)) warning("MARSSparamCIs: parSigma element has negative values on the diagonal.  These are replaced with NA.")
      
      vector.par.se[!is.neg] = try(sqrt(diag(MLEobj$parSigma)[!is.neg]), silent=TRUE)      # wrap in a try() in case it fails
      if(inherits(vector.par.se, "try-error") || any(is.nan(vector.par.se))) {
        warning("MARSSparamCIs: Some of the Hessian diagonals cannot be square-rooted; NA is being returned")
      }
    }
    vector.par.upCI = paramvec + qnorm(1-alpha/2)*vector.par.se
    vector.par.lowCI = paramvec - qnorm(1-alpha/2)*vector.par.se
    vector.par.bias = NULL; par.CI.nboot = NULL
  } #if method hessian
  
  if(method %in% c("parametric", "innovations"))  {
    boot.params = MARSSboot(MLEobj, nboot=nboot, output="parameters", sim=method,
                            param.gen="MLE", silent=silent)$boot.params
    vector.par.lowCI = apply(boot.params,1,quantile,probs=alpha/2)
    vector.par.upCI = apply(boot.params,1,quantile,probs=1-alpha/2)
    vector.par.se = sqrt(apply(boot.params,1,var))
    paramvec = MARSSvectorizeparam(MLEobj)
    vector.par.bias = paramvec - apply(boot.params,1,mean)
    par.CI.nboot = nboot
  }

  # Finalize the output
  names(vector.par.se)=paramnames
  if(!is.null(vector.par.bias)) names(vector.par.bias)=paramnames
  names(vector.par.upCI)=paramnames
  names(vector.par.lowCI)=paramnames
  for(val in c("par.se","par.bias","par.upCI","par.lowCI")){
    tmp.vec=get(paste("vector.",val,sep=""))
    if(!is.null(tmp.vec)) MLEobj[[val]] = MARSSvectorizeparam(MLEobj,tmp.vec)[["par"]]
    else MLEobj[[val]] = NULL
  }
  
  #add on information about how the CIs were constructed
  if(method=="hessian") MLEobj$par.CI.info = list(alpha=alpha, method=method, hessian.fun=hessian.fun) 
  if(method%in%c("parametric","innovations")) MLEobj$par.CI.info = list(alpha=alpha, method=method, nboot=par.CI.nboot) 
  return(MLEobj)
}

