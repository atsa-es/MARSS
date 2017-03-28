#######################################################################################################
#   MARSSmcinit function
#   Does a simple MonteCarlo initialization of the EM routine
#######################################################################################################
MARSSmcinit = function(MLEobj) {

  drawProgressBar = FALSE 
  if(!MLEobj$control$silent) { #then we can draw a progress bar
    cat("\n"); cat("> Starting Monte Carlo Initializations\n")
    prev = progressBar()
    drawProgressBar = TRUE
  }
  modelObj=MLEobj[["marss"]]
  y = modelObj$data
  par.dims=attr(modelObj,"model.dims")
  m = par.dims[["x"]][1]
  n = par.dims[["y"]][1]
  TT = par.dims[["data"]][2]
  ## YM matrix for handling missing values
  YM=matrix(as.numeric(!is.na(y)),n,TT)
  #Make sure the missing vals in y are zeroed out
  y[YM==0]=0

  free.tmp = modelObj$free 
  dim.tmp = list(Z=c(n,m), A=c(n,1), R=c(n,n), B=c(m,m), U=c(m,1), Q=c(m,m), x0=c(m,1))
  bounds.tmp = MLEobj$control$MCbounds
  init = bestinits = MLEobj$start
  bestLL = -1.0e10

  # loop over numInits: # of random draws of initial values
  for(loop in 1:MLEobj$control$numInits) {
    init.loop = init
      
    # Draw random values
    en = c("Z", "A", "R", "B", "U", "Q","x0")
    for(el in en) {
      dim.param = dim.tmp[[el]]
      if(!is.fixed(free.tmp[[el]])){
        bounds.param = bounds.tmp[[el]]
        #use the first fixed and free in a temporally varying model; arbitrary
        tmp=list(f=sub3D(modelObj$fixed[[el]],t=1),D=sub3D(modelObj$free[[el]],t=1))
        if(el %in% c("Q", "R")){   # random starts drawn from a wishart dist
	        if( bounds.param[1] < dim.param[1]){ df=dim.param[1] }else{ df=bounds.param[1] }
	        S=diag(bounds.param[2],dim.param[1])
	        #draw a random matrix from wishart
	        tmp.random = rwishart(df, S)/df
	        #reapply the sharing and fixed constraints 
          par.random = solve(t(tmp$D)%*%tmp$D)%*%t(tmp$D)%*%(vec(tmp.random)-tmp$f)
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
    MLEobj$control$maxit = MLEobj$control$numInitSteps
    MLEobj$control$minit = 1
    MLEobj$control$silent = TRUE #don't print convergence information during kem call          
    MLEobj = MARSSkem(MLEobj)  #get new fit using this init  

    if(drawProgressBar==TRUE) prev = progressBar(loop/MLEobj$control$numInits, prev)

    ## Check whether the likelihood is the best observed
## Revise: Only use bootstrap param draws where loglike did not go down during numInitSteps
    if(MLEobj$logLik > bestLL) {
      # update the best initial parameter estimates
      bestinits = MLEobj$par
      bestLL = MLEobj$logLik
    }
   
  } # end numInits loop

  return(bestinits)
}
