#######################################################################################################
#   MARSSparamCIs function
#   This returns CIs for ML parameter estimates
#######################################################################################################
MARSSparamCIs = function(MLEobj, method="hessian", alpha=0.05, nboot=1000, silent=TRUE) {
#this function expects a marssMLE object
#it will add some vectors of standard errors, biases, low/up CIs to MLEobj
if(!(method %in% c("hessian","parametric","innovations","hessian-harvey1989"))) stop("Stopped in MARSSparamCIs(). Current methods are hessian, parametric, innovations, and hessian-harvey1989.\n", call.=FALSE)
if(!is.marssMLE(MLEobj)) stop("Stopped in MARSSparamCIs(). This function needs a valid marssMLE object.\n", call.=FALSE)
paramvec = MARSSvectorizeparam(MLEobj)
if(length(paramvec)==0) stop("Stopped in MARSSparamCIs(). No estimated parameter elements.\n", call.=FALSE)
paramnames=names(paramvec)
if(method=="hessian")  {
    #Only for fixed, diagonal or block-diagonal unconstrained R, Q, or V0
    ok.form=TRUE
    RQV0type=c(Q=NA, R=NA, V0=NA)
    tmp=coef(MLEobj, type="matrix")
    for(elem in c("R","Q","V0")){
      if(length(MLEobj$par[[elem]])==0){ RQV0type[elem]="fixed"; next }  #ok is fixed
      if(is.diagonal(tmp[[elem]])){ RQV0type[elem]="diagonal"; next } #ok is diagonal
      if(grepl("unconstrained",MARSS:::describe.marssMODEL(MLEobj$marss)[[elem]])){
        RQV0type[elem]="unconstrained"; next 
        }  #ok is unconstrained
      ok.form==FALSE #hit an error
    }
    if(!ok.form) stop("MARSSparamCIs: To use method=hessian, R, Q, and V0 must be either fixed, diagonal or block diagonal unconstrained.")
    #Run emHessian to get hessian
    #MLEobj that is returned has var-cov matrices chol transformed
    MLEobj.hessian = MARSShessian(MLEobj)
    #deal with NAs in the Hessian
    na.diag = is.na(diag(MLEobj.hessian$Hessian))
    diag(MLEobj.hessian$Hessian)[ na.diag ] = 1 #this will lead to NA as needed
    MLEobj.hessian$Hessian[is.na(MLEobj.hessian$Hessian)]=0

    #standard errors
    hessInv = try(solve(MLEobj.hessian$Hessian))
    vector.par.se =rep(NA,length(paramvec)); names(vector.par.se)=names(paramvec)
    if(inherits(hessInv, "try-error") ) {
        warning("MARSSparamCIs: Hessian cannot be inverted; stderr = NA is being returned")
    }else{
      diag(hessInv)[ na.diag ] = 0
    is.neg=diag(hessInv)<0
    vector.par.se[!is.neg] = try(sqrt(diag(hessInv)[!is.neg]), silent=TRUE)      #invert the Hessian; wrap in a try() in case it fails
    if(inherits(vector.par.se, "try-error") || any(is.nan(vector.par.se))) {
        warning("MARSSparamCIs: Some of the hessian diagonals cannot be square-rooted; stderr = NA is being returned")
        }
    }    
    paramvec.hessian =  MARSSvectorizeparam(MLEobj.hessian)
    #this is in Cholesky transformed form
    vector.par.upCI = paramvec.hessian + qnorm(1-alpha/2)*vector.par.se
    vector.par.lowCI = paramvec.hessian - qnorm(1-alpha/2)*vector.par.se
    if(any(RQV0type=="diagonal")){ #need to back transform
      for(elem in c("Q","R","V0")[RQV0type=="diagonal"]){
        loc=grepl(elem, names(vector.par.upCI)) #where the parameter is in the vec
        tmp=MARSShessian.backtrans(MLEobj.hessian, vector.par.upCI)
        vector.par.upCI[loc] = tmp[loc]
        tmp=MARSShessian.backtrans(MLEobj.hessian, vector.par.lowCI)
        vector.par.lowCI[loc] = tmp[loc]
        #The inverse hessian we have is for the sqrt(var); we need the variance of the square of this
        #since the elem is diagonal, we can do the following
        var.elem=vector.par.se[loc]^2
        M1=MLEobj.hessian$parMean[loc] #first moment
        M2=M1^2+var.elem #second moment
        M4=M1^4+6*M1^2*var.elem+3*var.elem^2 #fourth moment
        vector.par.se[loc] = sqrt(M4-M2^2) #variance of elem^2 is M4-M2^2
      }
    }
    if(any(RQV0type=="unconstrained")){ #need to back transform with sim
      for(elem in c("Q","R","V0")[RQV0type=="unconstrained"]){
        loc=grepl(elem, names(vector.par.upCI)) #where the parameter is in the vec
        Mu=MLEobj.hessian$parMean
        Sigma=MLEobj.hessian$parSigma
        a=mvrnorm(1000,Mu,Sigma)
        elem.boot=matrix(NA,1000,sum(loc))
        for(i in 1:1000){
          tmp.MLE=MARSS:::MARSSvectorizeparam(MLEobj.hessian,a[i,])
          elem.chol=coef(tmp.MLE,type="matrix",form="marss")[[elem]]
          elem.tmp=elem.chol%*%t(elem.chol)
          elem.boot[i,]=elem.tmp[lower.tri(elem.tmp,diag=TRUE)]
        }
        tmp=apply(elem.boot,2,quantile,.975)
        vector.par.upCI[loc] = tmp
        tmp=apply(elem.boot,2,quantile,.025)
        vector.par.lowCI[loc] = tmp
        tmp=apply(elem.boot,2,var)
        vector.par.se[loc] = sqrt(tmp)
      }
    }
      vector.par.se[ na.diag ] = NA
      vector.par.lowCI[ na.diag ] = NA
      vector.par.upCI[ na.diag ] = NA
      vector.par.se[ na.diag ] = NA
      vector.par.bias = NULL; par.CI.nboot = NULL
} #if method hessian

if(method=="hessian-harvey1989")  {
  #Run Harvey 1989 recursion to get the hessian of the neg LL(observed FI matrix)
  MLEobj.hessian = MARSSFisherI(MLEobj)
  #deal with NAs in the Hessian
  na.diag = is.na(diag(MLEobj.hessian))
  diag(MLEobj.hessian)[ na.diag ] = 1 #this will lead to NA as needed
  MLEobj.hessian[is.na(MLEobj.hessian)]=0
  
  #standard errors
  hessInv = try(solve(MLEobj.hessian))
  vector.par.se =rep(NA,length(paramvec))
  if(inherits(hessInv, "try-error") ) {
    warning("MARSSparamCIs: Hessian cannot be inverted; stderr = NA is being returned")
  }else{
    diag(hessInv)[ na.diag ] = 0
    is.neg=diag(hessInv)<0
    vector.par.se[!is.neg] = try(sqrt(diag(hessInv)[!is.neg]), silent=TRUE)      #invert the Hessian; wrap in a try() in case it fails
    if(inherits(vector.par.se, "try-error") || any(is.nan(vector.par.se))) {
      warning("MARSSparamCIs: Some of the hessian diagonals cannot be square-rooted; stderr = NA is being returned")
    }
  }
  #Harvey recursion is a hessian for non-transformed parameters
  paramvec.hessian = paramvec
  vector.par.upCI = paramvec.hessian + qnorm(1-alpha/2)*vector.par.se
  vector.par.lowCI = paramvec.hessian - qnorm(1-alpha/2)*vector.par.se
  vector.par.lowCI[ na.diag ] = NA
  vector.par.upCI[ na.diag ] = NA
  vector.par.se[ na.diag ] = NA
  vector.par.bias = NULL; par.CI.nboot = NULL
} #if method hessian-harvey1989

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
names(vector.par.se)=paramnames
if(!is.null(vector.par.bias)) names(vector.par.bias)=paramnames
names(vector.par.upCI)=paramnames
names(vector.par.lowCI)=paramnames
for(val in c("par.se","par.bias","par.upCI","par.lowCI")){
  tmp.vec=get(paste("vector.",val,sep=""))
  if(!is.null(tmp.vec)) MLEobj[[val]] = MARSSvectorizeparam(MLEobj,tmp.vec)[[val]]
  else MLEobj[[val]] = NULL
}
MLEobj$par.se = MARSSvectorizeparam(MLEobj,vector.par.se)[["par"]]
if(!is.null(vector.par.bias)) MLEobj$par.bias = MARSSvectorizeparam(MLEobj,vector.par.bias)[["par"]]
  else MLEobj$par.bias = NULL
MLEobj$par.upCI = MARSSvectorizeparam(MLEobj,vector.par.upCI)[["par"]]
MLEobj$par.lowCI = MARSSvectorizeparam(MLEobj,vector.par.lowCI)[["par"]]
MLEobj$par.CI.alpha = alpha
MLEobj$par.CI.method = method
MLEobj$par.CI.nboot = par.CI.nboot
return(MLEobj)
}

