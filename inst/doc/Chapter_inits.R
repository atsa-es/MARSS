###################################################
### code chunk number 2: Cs_mci_001
###################################################
fulldat = lakeWAplanktonTrans
years = fulldat[,"Year"]>=1965 & fulldat[,"Year"]<1975
dat = t(fulldat[years,c("Greens", "Bluegreens")])
the.mean = apply(dat,1,mean,na.rm=TRUE)
the.sigma = sqrt(apply(dat,1,var,na.rm=TRUE))
dat = (dat-the.mean)*(1/the.sigma)


###################################################
### code chunk number 3: Cs_mci_002
###################################################
covariates = rbind(
   Temp = fulldat[years,"Temp"],
   TP = fulldat[years,"TP"])
# demean the covariates
the.mean = apply(covariates,1,mean,na.rm=TRUE)
covariates = covariates-the.mean


###################################################
### code chunk number 4: Cs_mci_003
###################################################
U = x0 = "zero"
Q = "unconstrained"
d = covariates
A = "zero"
D = "unconstrained"
R = "diagonal and equal"
model.list = list(U=U,Q=Q,A=A,R=R,D=D,d=d,x0=x0)
kem = MARSS(dat, model=model.list)


###################################################
### code chunk number 5: Cs_mci_004
###################################################
coef(kem, what="par")


###################################################
### code chunk number 6: Cs_mci_005
###################################################
out=coef(kem, what="par")
out$D
out$Q


###################################################
### code chunk number 7: Cs_mci_006
###################################################
inits=list(Q=1)
kem = MARSS(dat, model=model.list, inits=inits)
#or
inits=list(Q=matrix(c(1,0,1),3,1))
kem = MARSS(dat, model=model.list, inits=inits)


###################################################
### code chunk number 8: Cs_mci_007
###################################################
inits=list(Q=matrix(c(1,0.5,0.7),3,1))
kem = MARSS(dat, model=model.list, inits=inits)


###################################################
### code chunk number 9: Cs_mci_008
###################################################
inits=list(Q=matrix(c(1,0.5,0.7),3,1), D=1)
kem = MARSS(dat, model=model.list, inits=inits)


###################################################
### code chunk number 10: Cs_mci_009
###################################################
inits=list(D=coef(kem, what="par")$D)
kem = MARSS(dat, model=model.list, inits=inits)


###################################################
### code chunk number 11: Cs_mci_010
###################################################
#create the par list from the output
inits=coef(kem, what="par")
kem.bfgs = MARSS(dat, model=model.list, inits=inits, method="BFGS")


###################################################
### code chunk number 12: Cs_mci_0101
###################################################
#create the par list from the output
kem.bfgs = MARSS(dat, model=model.list, inits=kem, method="BFGS")


###################################################
### code chunk number 13: Cs_mci_011
###################################################
######################################################################################################   MARSSmcinit function
#   Does a simple MonteCarlo initialization of the EM routine
#   The function uses a number of MARSS utility functions accessed with MARSS:::
#######################################################################################################
MARSSmcinit = function(MLEobj, 
                       control=list(numInits = 500, numInitSteps = 10,
                                    MCbounds=list(B=c(0,1), U=c(-1,1), Q = c(1,1), 
                                                  Z=c(0,1), A=c(-1,1), R = c(1,1), x0 = c(-1,1) )),
                       silent=FALSE) {
  
  control.default=list(numInits = 500, numInitSteps = 10, MCbounds=list(B=c(0,1), U=c(-1,1), Q = c(1,1), Z=c(0,1), A=c(-1,1), R = c(1,1), x0 = c(-1,1) ))
  if(!is.null(control)){
    if(!is.list(control)) stop("MARSSmcinit: control must be a list")
    if(any(!(names(control)%in%names(control.default)))) stop(paste("MARSSmcinit: allowed control list elements are", names(control.default)))
    control.new=control.default
    for(i in names(control)) control.new[[i]]=control[[i]]
    control = control.new
  }
  drawProgressBar = FALSE 
  if(!silent) { #then we can draw a progress bar
    cat("\n"); cat("> Starting Monte Carlo Initializations\n")
    prev = MARSS:::progressBar() #this is an internal function to MARSS
    drawProgressBar = TRUE
  }
  MODELobj=MLEobj[["marss"]]
  y = MODELobj$data
  par.dims=attr(MODELobj,"model.dims")
  m = par.dims[["x"]][1]
  n = par.dims[["y"]][1]
  TT = par.dims[["data"]][2]
  ## YM matrix for handling missing values
  YM=matrix(as.numeric(!is.na(y)),n,TT)
  #Make sure the missing vals in y are zeroed out
  y[YM==0]=0

  free.tmp = MODELobj$free 
  dim.tmp = list(Z=c(n,m), A=c(n,1), R=c(n,n), B=c(m,m), U=c(m,1), Q=c(m,m), x0=c(m,1))
  bounds.tmp = control$MCbounds
  init = bestinits = MLEobj$start
  bestLL = -1.0e10

  # loop over numInits: # of random draws of initial values
  for(loop in 1:control$numInits) {
    init.loop = init
      
    # Draw random values
    en = c("Z", "A", "R", "B", "U", "Q","x0")
    for(el in en) {
      dim.param = dim.tmp[[el]]
      if(!MARSS:::is.fixed(free.tmp[[el]])){ #is.fixed is a utility func in MARSS
        bounds.param = bounds.tmp[[el]]
        #use the first fixed and free in a temporally varying model; arbitrary
        tmp=list(f=MARSS:::sub3D(MODELobj$fixed[[el]],t=1),D=MARSS:::sub3D(MODELobj$free[[el]],t=1))
        if(el %in% c("Q", "R")){   # random starts drawn from a wishart dist
          if( bounds.param[1] < dim.param[1]){ df=dim.param[1] }else{ df=bounds.param[1] }
	        S=diag(bounds.param[2],dim.param[1])
	        #draw a random matrix from wishart
	        tmp.random = MARSS:::rwishart(df, S)/df
	        #reapply the sharing and fixed constraints 
          par.random = solve(t(tmp$D)%*%tmp$D)%*%t(tmp$D)%*%(MARSS:::vec(tmp.random)-tmp$f)
	      }else{
	        par.random = matrix(runif(dim(tmp$D)[2], bounds.param[1], bounds.param[2]), dim(tmp$D)[2],1)
          if(el %in% c("B")){
           tmp.max=max(abs(eigen(par.random,only.values=TRUE)$values))
           #rescale to bring the max abs eigenvalues to between 0 and 1
           par.random =  par.random/(tmp.max/runif(1,.01,.99))
          }
          if(el %in% c("x0")){
           x0init = init$x0 #where the original start is
           x.lo = ifelse(x0init > 0, exp(bounds.param[1])*x0init, exp(bounds.param[2])*x0init)
           x.hi = ifelse(x0init > 0, exp(bounds.param[2])*x0init, exp(bounds.param[1])*x0init)
           par.random = matrix(runif(dim(tmp$D)[2], x.lo, x.hi), dim(tmp$D)[2],1)
          }

        } 
      }else{ par.random=matrix(0,0,1) }
      init.loop[[el]] = par.random 
    }

    ## Call MARSSkem() with these inits 
    MLEobj$start = init.loop
    MLEobj$control$maxit = control$numInitSteps
    MLEobj$control$minit = 1
    MLEobj$control$silent = TRUE #don't print convergence information during kem call          
    MLEobj = MARSSkem(MLEobj)  #get new fit using this init  

    if(drawProgressBar==TRUE) prev = MARSS:::progressBar(loop/control$numInits, prev)

    ## Check whether the likelihood is the best observed
    ## Only use bootstrap param draws where loglike did not go down during numInitSteps
    if(MLEobj$logLik > bestLL) {
      # update the best initial parameter estimates
      bestinits = MLEobj$par
      bestLL = MLEobj$logLik
    }
   
  } # end numInits loop

  return(bestinits)
}


###################################################
### code chunk number 14: Cs_mci_012
###################################################
  dat = t(harborSeal)
  dat = dat[c(2,nrow(dat)),]
  fit1=MARSS(dat)
  MCinits = MARSSmcinit(fit1, control=list(numInits = 10)) 
  fit2=MARSS(dat, inits=MCinits)


