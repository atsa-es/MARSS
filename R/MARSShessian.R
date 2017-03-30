#######################################################################################################
#   MARSShessian.chol functions
#   Numerical estimation of the hessian of the negative log-likelihood function
#   Uses a cholesky transformation of the var-cov matrices to ensure they stay positive definite
#
#   Adds Hessian, parameter var-cov matrix, and parameter mean to a marssMLE object
#   Returns par and fixed and free in chol form
#   The func MARSShessian.chol.backtrans (below) transforms that to the original form
#######################################################################################################
MARSShessian.chol = function(MLEobj) {
  ##
  # Define needed functions
  ##
  kfNLL = function(x, MLEobj=NULL){  #NULL assignment needed for optim call syntax
    #MLEobj is the MLEobj with chol transformation of the variances so has altered free and fixed
    #x is the paramvector with variances chol transformed
    
    par.untransformed=MARSShessian.backtrans(MLEobj, x)
    MLEobj.orig=MLEobj
    MLEobj.orig$marss$fixed = MLEobj$fixed.original
    MLEobj.orig$marss$free = MLEobj$free.original
    #reconstruction the original
    MLEobj.orig=MARSSvectorizeparam(MLEobj.orig, par.untransformed)
    
    #kfsel selects the Kalman filter / smoother function based on MLEobj$fun.kf
    negLL = MARSSkf(MLEobj.orig, only.logLik=TRUE, return.lag.one=FALSE )$logLik

    -1*negLL
  }
  
  
  fun="MARSSkf"
  
  ## attach would be risky here since user might have one of these variables in their workspace    
  y = MLEobj[["marss"]][["data"]] #must have time going across columns
  marss.object = MLEobj[["marss"]]
  free = marss.object[["free"]]
  fixed = marss.object[["fixed"]]
  par.dims=attr(marss.object,"model.dims")
  pars=MLEobj[["par"]]
  
  #The code is used to set things up to use MARSSvectorizeparam to just select inits for the estimated parameters
  tmp.MLEobj = MLEobj
  #This is needed for the likelihood calculation
  tmp.MLEobj$fixed.original = tmp.MLEobj$marss$fixed
  tmp.MLEobj$free.original = tmp.MLEobj$marss$free
  
  for(elem in c("Q","R","V0")){ #need the chol for these
    if(!is.fixed(free[[elem]])){
      
      tmp.par=matrix(0,dim(tmp.MLEobj[["par"]][[elem]])[1],1) #holder for the estimated elements
      TT.f=dim(fixed[[elem]])[3]
      TT.d=dim(free[[elem]])[3]
      for(t in 1:max(TT.f,TT.d)){
        #This requires a chol transformation and that trick only works for certain var-cov structures.
        #method="BFGS" has the same constraints so I can use is.validvarcov() to test
        #test each var-cov matrix at each time step in case time-varying
        par.as.list = fixed.free.to.formula(fixed[[elem]][,,min(TT.f,t),drop=FALSE],free[[elem]][,,min(TT.d,t),drop=FALSE], par.dims[[elem]][1:2]) #coverts the fixed,free pair to a list matrix
        tmp=is.validvarcov(par.as.list, method="BFGS")
        if(!tmp$ok) stop("Stopped in MARSShessian(): The variance matrix must be diagonal for the Hessian computation.")#I think you can have time-varying but I need to figure how to compute fixed and free using the t=1 par; 
        
        #what's happening here is I am chol-transforming the var-cov matrix.  The matrix can have fixed and shared elements.
        #This is why I can't just do chol(the.par).  I need to solve for the pars
        f=sub3D(fixed[[elem]],t=min(t,TT.f)) 
        d=sub3D(free[[elem]],t=min(t,TT.d))
        the.par=unvec(f+d%*%pars[[elem]], dim=par.dims[[elem]][1:2])
        is.zero=diag(the.par)==0   #where the 0s on diagonal are
        if(any(is.zero)) diag(the.par)[is.zero]=1    #so the chol doesn't fail if there are zeros on the diagonal
        the.par=t(chol(the.par))  #transpose of chol
        if(any(is.zero)) diag(the.par)[is.zero]=0  #set back to 0
        
        #when being passed to fdHess, pars for var-cov mat is the chol which has the upper.tri set to 0, 
        #so need to reset free and fixed matrices
        #compute the D matrix corresponding to upper.tri=0 at in t(chol)
        #This part that only works if (block) diagonal or (block) unconstrained; in only those cases the 
        #is list matrix for the chol simply the list matrix for original with 0s in upper tri
        tmp.list.mat=fixed.free.to.formula(sub3D(tmp.MLEobj[["marss"]][["fixed"]][[elem]],t=min(t,TT.f)),sub3D(tmp.MLEobj[["marss"]][["free"]][[elem]],t=min(t,TT.d)),par.dims[[elem]][1:2])
        tmp.list.mat[upper.tri(tmp.list.mat)]=0   #set upper tri to zero
        d.tchol=convert.model.mat(tmp.list.mat)[["free"]][,,1] #because only sending one f and d to fixed.free.to.formulat
        tmp.MLEobj[["marss"]][["free"]][[elem]][,,min(t,TT.d)]=d.tchol #set the correct d[,,t] to the chol version
        #I don't need the fixed matrix at all
        
        #This part also doesn't work if not (block) diagonal or (block) unconstrained
        #This only works because I don't allow matrix elements to have 2 estimated values or estimated and fixed values; a+b and 1+a are illegal
        #the diag(as.numeric(tmp.par==0)) is removing par estimates that have already been computed (!=0) in prev time steps
        #from f+Dm=M so m = solve(crossprod(d))%*%t(d)%*%(vec(the.par)-f)
        #but if d!=0,then f==0. if f!-0, then d==0.  
        #Thus crossprod(d))%*%t(d) has 0 cols where fs appear in the.par and f is not needed
        tmp.par = tmp.par + diag(as.numeric(tmp.par==0))%*%solve(crossprod(d.tchol))%*%t(d.tchol)%*%vec(the.par)
      }
      tmp.MLEobj[["par"]][[elem]] = tmp.par 

    }else{ tmp.MLEobj[["par"]][[elem]] = matrix(0,0,1) }

  }
  # will return the inits only for the estimated parameters
  paramvector = MARSSvectorizeparam(tmp.MLEobj)
  
  #this is the one that is chol transformed
  MLEobj=tmp.MLEobj
  
  #   kfNLL=function(paramvec, MLEobj){
  #     #neglogLik is defined in MARSSoptim
  #     return(neglogLik(paramvec, MLEobj=MLEobj))
  #   }
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
  
  #This is the TRANSFORMED MLEobj
  #also has fixed.original and free.original
  return(MLEobj)
}

MARSShessian.chol.backtrans = function(MLEobj.hessian, par.hessian){
  #MLEobj is your original untransformed MLEobj
  #MLEobj.hessian is a transformed version with the chol transformation for variances
  #par.hessian is a vector of parameters where the variances are in chol form
  #Goal is to put the par variances values in 
  
  #first put the parameters into the MLEobj.hession object
  MLEobj.hessian=MARSSvectorizeparam(MLEobj.hessian, par.hessian)
  free = MLEobj.hessian[["marss"]][["free"]]
  pars = MLEobj.hessian[["par"]]
  
  #chol back transformation
  par.dims=attr(MLEobj.hessian[["marss"]],"model.dims")
  for(elem in c("Q","R","V0")){   #this works because by def fixed and free blocks of var-cov mats are independent
    if(!is.fixed(free[[elem]])) #if not estimated then there won't be a par element
    {
      tmp.par=matrix(0,dim(pars[[elem]])[1],1) #holder for the estimated elements
      TT.d=dim(free[[elem]])[3]
      for(t in 1:TT.d){
        d=sub3D(free[[elem]],t=t)
        par.dim=par.dims[[elem]][1:2]
        #t=1 since D not allowed to be time-varying; since code 4 lines down won't work otherwise
        L=unvec(d%*%pars[[elem]],dim=par.dim) #this by def will have 0 row/col at the fixed values
        the.par = tcrossprod(L) #L%*%t(L)
        tmp.par = tmp.par + diag(as.numeric(tmp.par==0))%*%solve(crossprod(d))%*%t(d)%*%vec(the.par)
      }
      MLEobj.hessian[["par"]][[elem]] = tmp.par 
    }
    
  } #end for
  
  #now the MLEobj par elements are back transformed
  return( MARSSvectorizeparam(MLEobj.hessian) )
}
