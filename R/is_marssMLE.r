#####
# is.marssMLE()
# Check that the marssMLE object has all the parts it needs 
# and that these have the proper size and form
#####
is.marssMLE <- function(MLEobj) 
{
  if( !("marssMLE" %in% class(MLEobj)) ) stop("Stopped in is.marssMLE() because object class is not marssMLE.\n", call.=FALSE)
  
  msg = c()
 ## Check for required components
 el = c("marss","model", "start", "control", "method")
 if( !all(el %in% names(MLEobj)) ){
    msg = c(msg, paste("Element", el[!(el %in% names(MLEobj))], "is missing from object.\n"))
 }
  ## Break out now if there was a problem
  if(length(msg) != 0){ 
    msg=c("\nErrors were caught in is.marssMLE()\n", msg)
    return(msg)
  }
  
  #here form="marss" since $marss is of that form
  msg = is.marssMODEL(MLEobj[["marss"]], method=MLEobj[["method"]])        #returns TRUE or a vector of msgs
  ## Break out now if there was a problem with the model
  if(!isTRUE(msg)){ 
    msg=c("\nErrors were caught in is.marssMLE() in the call to is.marssMODEL().\n", msg)
    return(msg)
  }

  #check that model object in the called form is ok too
  msg = is.marssMODEL(MLEobj[["model"]], method=MLEobj[["method"]])        #returns TRUE or a vector of msgs
  ## Break out now if there was a problem with the model
  if(!isTRUE(msg)){ 
    msg=c("\nErrors were caught in is.marssMLE() in the call to is.marssMODEL().\n", msg)
    return(msg)
  }
  
  ## is.wholenumber() borrowed from is.integer example
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

  msg = c()   #reset the messages

  ## If model is OK, check start
  model=MLEobj[["marss"]]
  free=model[["free"]]
  fixed=model[["fixed"]]
  par.dims=attr(model,"model.dims")
  en=attr(model,"par.names")
  dat = model[["data"]]

  init.null = dim.init = NULL

  #make sure each element of start is present and is numeric
  for (el in en) {   
    init.null.flag <- ( is.null(MLEobj[["start"]][[el]]) || !is.numeric(MLEobj[["start"]][[el]]) )

    dim.init.flag = FALSE 
    if (!init.null.flag) {   #element is present so check it's dimensions
      #true means there is a problem
      dim.init.flag = (dim(free[[el]])[2]!=dim(MLEobj[["start"]][[el]])[1]) || dim(MLEobj[["start"]][[el]])[2]!=1
    }
 
    init.null <- c(init.null, init.null.flag)
    dim.init <- c(dim.init, dim.init.flag)  
  }
    problem <- any(c(init.null, dim.init))
  if (problem) {    
    if(any(init.null)) {
      msg = c(msg, paste("Missing or non-numeric initial value", en[init.null],"\n"))
    }
    if(any(dim.init)) {
      msg = c(msg, paste("Dims for initial value ", en[dim.init], " do not match model specification.\n", sep=""))
    }    
  }
## TEMPORARY UNTIL change Q and R to not allow 0s (with addition of G and H matrices)
  #make sure not 0s in par unless corresponding rows of D are all zero
  for(el in c("Q","R","V0")){
  if(!is.fixed(free[[el]])){
    bad.ts=c()
    dvars=free[[el]][1 + 0:(par.dims[[el]][1] - 1)*(par.dims[[el]][1] + 1),,,drop=FALSE]
    dvars.cols=apply(dvars!=0,2,sum)!=0   #those var rows that have values, give use the var columns
    Tmax=max(dim(fixed[[el]])[3],dim(free[[el]])[3])
    for(i in 1:Tmax){
      d=sub3D(free[[el]],t=min(i,dim(free[[el]])[3]))
      if(any(colSums(d[,dvars.cols,drop=FALSE])!=0 & MLEobj[["start"]][[el]][dvars.cols]==0)){
        bad.ts=c(bad.ts,i)
      }
    }
    if(length(bad.ts)!=0){
     if(length(bad.ts)<5){ ts.char=paste(bad.ts,collapse=", ") 
     }else{ ts.char=paste(paste(bad.ts[1:5], collapse=", "),",...") }
     msg = c(msg, paste("At t=",ts.char,", one of the start$", el, " for a variance is 0,\n and the corresponding column of model$free$",el,"[,,t] matrix is not all zero.\n",sep=""))
     }
  }
  }

  ## Check params consistency if present
  if(!is.null(MLEobj[["par"]])) {
    dim.par = NULL
    for (el in en) {
      dim.flag=FALSE
      if(!is.null(MLEobj[["par"]][[el]])) {
         dim.flag = isTRUE(dim(MLEobj[["par"]][[el]])[2]!=1)
         dim.flag = dim.flag || isTRUE(dim(MLEobj[["par"]][[el]])[1]!=dim(free[[el]])[2])
         dim.par=c(dim.par,dim.flag)
      }
    }
    if(any(dim.par)){
       msg = c(msg, paste("par element for", en[dim.par],"is incorrect.  Should correspond to dim 2 of free.\n"))
    }
  }
  
  alldefaults=get("alldefaults", envir=pkg_globals)
  ## Check controls
  if(!is.null(MLEobj[["control"]])){
    if(!is.list(MLEobj[["control"]])) stop("Stopped in is.marssMLE() because control must be passed in as a list.\n", call.=FALSE)
  control = MLEobj[["control"]]
  en = names(alldefaults[[MLEobj[["method"]]]][["control"]])
  ok.null=c("REPORT", "reltol", "fnscale", "parscale", "ndeps", "alpha", "beta", "gamma", "type", "lmm", "factr",
      "pgtol", "tmax", "temp", "lower", "upper")
  for (el in en) {
    null.flag <- ( is.null(control[[el]]) && !(el %in% ok.null ) )  #those in ok.null can be NULL
    if(null.flag) msg = c(msg, paste(el,"is missing from the control list\n"))

    if( !is.null(control[[el]]) ) { 
      #everything must be numeric except these
      if( el %in% en[!(en %in% c("safe", "allow.degen", "demean.states", "sparse"))] ) {
        null.flag <- (!is.numeric(control[[el]]))
        if(null.flag) msg = c(msg, paste("control list element", el,"is non-numeric\n"))
      }
      #everything must be positive or 0 except these
      if( (el %in% en[!(en %in% c("safe", "trace", "allow.degen", "demean.states", "sparse"))])  && is.numeric(control[[el]]) ) {
        null.flag <- ( control[[el]] <= 0)
        if(null.flag) msg = c(msg, paste("control list element", el,"less than or equal to zero\n"))
      }
      #trace can be -1 but no smaller
      if ( (el=="trace") && is.numeric(control[[el]]) ) {
        null.flag <- ( control[[el]] < -1)
        if(null.flag) msg = c(msg, paste("control list element trace must be an integer greater than or equal to -1.\n"))
      }
      #these need to be whole numbers
      if (el %in% c("trace", "minit", "maxit", "min.iter.conv.test", "conv.test.deltaT", "min.degen.iter") && is.numeric(control[[el]])) {
        null.flag <- ( !is.wholenumber(control[[el]]) )
        if(null.flag) msg = c(msg, paste("control list element", el,"is not a whole number\n"))
      }
      #minit must be less than maxit
      if (el %in% c("minit") && !is.null(control[["maxit"]]) ) {
        if(!is.null(control[["minit"]]) && !is.null(control[["maxit"]]) && is.numeric(control$minit) && is.numeric(control$maxit) && is.wholenumber(control$minit) &&  is.wholenumber(control$maxit)) null.flag <- (control$minit > control$maxit)  
        if(null.flag) msg = c(msg, paste("control list element minit is greater than maxit\n"))
      }
      #conv.test.deltaT can't be less than 2
      if (el %in% c("conv.test.deltaT") && is.numeric(control[[el]]) ) {
        null.flag <- ( control[[el]] < 2)
        if(null.flag) msg = c(msg, "control list element conv.test.deltaT must be greater than 2\n")
      }
      if (el %in% c("safe", "allow.degen", "demean.states", "sparse")) {
        null.flag <- !(control[[el]] %in% c(TRUE, FALSE) )	  
        if(null.flag) msg = c(msg, paste("control list element", el,"is not TRUE or FALSE\n"))
      }

    } # el is not null     
  }      # for el in en
  } #not null control
  
  ## if $diffuse, method="BFGS" and $tinitx=1
  if( identical(MLEobj$marss$diffuse,TRUE) & !(MLEobj$method=="BFGS" & MLEobj$marss$tinitx==1) ) 
    msg = c(msg, "If you specify a diffuse prior, method must be BFGS and model$tinitx set to 1.\n")    
 

if(length(msg) == 0){ return(TRUE)
}else {
  msg=c("\nErrors were caught in is.marssMLE(). Type MARSSinfo() for more information.\n", msg)
  return(msg)
}
}
