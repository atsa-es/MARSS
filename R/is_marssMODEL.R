########################################################################
# is.marssMODEL function
# Check that the MODELobj object has all the parts it needs
# data, fixed, free, and X.names
# and that these have the proper size and form
# That it has all its attributes and 
# that the fixed and free matrices fit the specified attributes (internally consistent)
########################################################################
is.marssMODEL <- function(MODELobj, method="kem") {
  if( !inherits(MODELobj ,"marssMODEL") ) stop("Stopped in is.marssMODEL() because object class is not marssMODEL.\n", call.=FALSE)
  msg = NULL
  
  ###########################
  # First do some basic mode and presence tests so that the other tests work
  ###########################
  
  ## Check for required components
  el = c("data","fixed","free","tinitx","diffuse")
  if( !all(el %in% names(MODELobj)) ) { 
    msg = c(msg, "Element", el[!(el %in% names(MODELobj))], "is missing from the model object.\n")
  }
  if(!is.null(msg)){  #rest of the tests won't work so stop now
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }
  
  ## Check that free and fixed are numeric matrices with no NA or Infs
  for(mat in c("fixed","free")){
    if (!is.list(MODELobj[[mat]])) msg = c(msg, paste("MODELobj$",mat," must be a list of matrices.\n",sep="")) 
    for (i in 1:length(MODELobj[[mat]])) {
      if( !is.array(MODELobj[[mat]][[i]]) && length(dim(MODELobj[[mat]][[i]]))!=3 ){ 
        msg = c(msg, paste("MODELobj$",mat,"$",names(MODELobj[[mat]])[i]," must be a 3D matrix.\n", sep=""))
      }
      if(mode(MODELobj[[mat]][[i]]) != "numeric" || any(is.na(MODELobj[[mat]][[i]])) || any(is.infinite(MODELobj[[mat]][[i]])) ) 
        msg = c(msg, paste("MODELobj$",mat,"$",names(MODELobj[[mat]])[i]," must be numeric, have no NAs, and no Infs.\n", sep=""))  
    }
  }
  
  if( length(dim(MODELobj$data)) != 2)
    msg = c(msg, "Data is not a 2D matrix.\n")
  ## check for T=1
  if( !is.numeric(MODELobj$data ) ) msg = c(msg, "Data must be numeric.\n")
  if( dim(MODELobj$data)[2] == 1 ) msg = c(msg, "Data has only one time point.\n")
  
  if(!is.null(msg)){  #rest of the tests won't work so stop now
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }
  
  ###########################
  # Check that the attributes are complete and consistent
  ###########################
  el = c("form","model.dims","par.names","X.names","Y.names", "equation","obj.elements")
  attr.names=names(attributes(MODELobj))
  if( !all(el %in% attr.names) ) { 
    msg = c(msg, "Element", el[!(el %in% attr.names)], "is missing from the attributes of the model object.\n")
  }
  if(!is.null(msg)){  #rest of the tests won't work so stop now
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }  
  el=c("par.names","form","X.names","Y.names","equation")
  for(elem in el){
    fattr=attr(MODELobj,elem)
    if(!is.vector(fattr) || !is.character(fattr)){
      msg=c("The ", elem, " attribute of the marssMODEL object needs to be a character vector.\nErrors were caught in is.marssMODEL()\n", msg, sep="")
      return(msg)   #the rest of the tests will hang so stop now 
    }
  }
  
  par.names=attr(MODELobj,"par.names")
  if(any(duplicated(par.names)))
    msg=c(msg, "par.names attribute of the model object has duplicated names.\n")
  ###########################
  # Check that fixed and free have all names in par.names
  # and that no names in fixed / free that aren't in par.names
  ###########################
  el = c("fixed","free")
  for(elem in el){
    fnames=names(MODELobj[[elem]])
    if( !all(par.names %in% fnames ) ) { 
      msg = c(msg, "Element ", par.names[!(par.names %in% fnames)], " is missing from the ", elem," element of the model object.\n", sep="")
    }
    if( !all( fnames %in% par.names ) ) { 
      msg = c(msg, "Element ", fnames[!(fnames %in% par.names)], " in ", elem, " is missing from the par.names attribute of the model object.\n", sep="")
    }
  }

  ###########################
  # Check that model dims have all the par.names
  ###########################
  model.dims=attr(MODELobj,"model.dims")
  #the info in model.dims should be par.names and the extras
  model.dim.names=c(par.names,"data","x","y","w","v")
  fnames=names(model.dims)
  fnames=fnames[!(fnames %in% c("data","x","y","w","v"))]
  if(!is.list(model.dims))
    msg=c(msg,"model.dims attribute must be a list.\n")
  if( !all(par.names %in% fnames) ) { 
    msg = c(msg, "Element", par.names[!(par.names %in% fnames)], "is missing from the model.dims attribute of the model object.\n", sep="")
  } 
  if( !all( fnames %in% par.names ) ) { 
    msg = c(msg, "Element ", fnames[!(fnames %in% par.names)], " in model.dims attribute is missing from the par.names attribute of the model object.\n", sep="")
  }
  if(!is.null(msg)){  #rest of the tests won't work so stop now
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }
  ###########################
  # Check that length of X.names matches first dim of model.dims$X
  # Check that length of Y.names matches first dim of model.dims$Y
  ###########################
  if(length(attr(MODELobj,"X.names"))!=attr(MODELobj,"model.dims")$x[1])
    msg="The length of the X.names attribute of model object must equal the first element of the model.dims attribute x element.\n"
  if(length(attr(MODELobj,"Y.names"))!=attr(MODELobj,"model.dims")$y[1])
    msg=c(msg,"The length of the Y.names attribute of model object must equal the first element of the model.dims attribute y element.\n")
  ###########################
  # Check that 2nd dim of model.dims$x and model.dims$y equals the 2nd dim of data
  ###########################
  TT=dim(MODELobj$data)[2]
  if(attr(MODELobj,"model.dims")$y[2]!=TT || attr(MODELobj,"model.dims")$x[2]!=TT ||
       attr(MODELobj,"model.dims")$w[2]!=TT || attr(MODELobj,"model.dims")$v[2]!=TT )
    msg=c(msg,"The 2nd element of the model.dims attribute for x, y, w, and v must equal the number of time points in the data.\n")

  if(!is.null(msg)){  #rest of the tests won't work so stop now
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }
  
  ###########################
  # Check that fixed and free are 3D
  ###########################
  dim.fixed = dim.free =  NULL
  free=MODELobj$free; fixed=MODELobj$fixed
  for (elem in par.names) {
    dim.fixed.flag = dim.free.flag = FALSE
    if(length(dim(free[[elem]]))!=3){ dim.free.flag = TRUE } #3-dimensions
    if(length(dim(fixed[[elem]]))!=3){ dim.fixed.flag = TRUE }
    dim.fixed <- c( dim.fixed, dim.fixed.flag )
    dim.free <- c( dim.free, dim.free.flag )
  }
  #stop now since the rest of the tests won't work
  if (any(c(dim.fixed, dim.free))) {  
    if(any(dim.fixed)) {
      msg = c(msg, paste("fixed", par.names[dim.fixed],"is not 3D.\n"))   }
    if(any(dim.free)) {
      msg = c(msg, paste("free", par.names[dim.free],"is not 3D.\n"))    }
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }

  ###########################
  # Check that 1st and 2nd dims of fixed ok, and 1st dim of free are ok
  ###########################
  dim.fixed.flag = dim.free.flag = nomatch.flag = FALSE
  dim.fixed = dim.free = NULL
  
  for (elem in par.names) {
    ## Check for problems in the fixed/free pairs. Problems show up as TRUE     
    # check dim
    dim.fixed.flag = !isTRUE(all.equal( dim(fixed[[elem]])[1], model.dims[[elem]][1]*model.dims[[elem]][2] ) )
    dim.fixed.flag = dim.fixed.flag || !isTRUE( all.equal( dim(fixed[[elem]])[2], 1 ) )    
    dim.free.flag = !isTRUE(all.equal( dim(free[[elem]])[1], model.dims[[elem]][1]*model.dims[[elem]][2] ) )
    dim.fixed = c(dim.fixed, dim.fixed.flag)
    dim.free = c(dim.free, dim.free.flag) 
  }
  if (any(c(dim.fixed, dim.free))) {  #There's a problem
    if(any(dim.fixed)) {
      bad.names = par.names[dim.fixed]
      msg = c(msg, paste("fixed ", bad.names, " dims are incorrect. Dims 1 and 2 should be (", unlist(lapply(model.dims[bad.names],function(x){x[1]})), "x", unlist(lapply(model.dims[bad.names],function(x){x[2]})),", 1) based on data and other parameters.\n",sep=""))
    }
    if(any(dim.free)) {
      bad.names = par.names[dim.free]
      msg = c(msg, paste("free", bad.names, "dims are incorrect. Dim 1 be ", model.dims[bad.names], "x", model.dims[bad.names],"based on data and other parameters.\n"))
    }
  }

  ###########################
  # Check that 3rd dims of fixed and free are 1 or TT
  ###########################
  dim.fixed=dim.free=c()
  dim.fixed.flag = dim.free.flag  = FALSE
  for (elem in par.names) {
    ## Check for problems in the fixed/free pairs. Problems show up as TRUE 
    
    #test that dim3 is either 1 or TT
    dim.fixed.flag = (!isTRUE( all.equal( dim(fixed[[elem]])[3], TT ) ) & !isTRUE(all.equal( dim(fixed[[elem]])[3], 1 ) )) 
    dim.free.flag = (!isTRUE(all.equal( dim(free[[elem]])[3], TT ) ) & !isTRUE(all.equal( dim(free[[elem]])[3], 1 ) ))    
    dim.fixed = c(dim.fixed, dim.fixed.flag)
    dim.free = c(dim.free, dim.free.flag) 
  }  

  if (any(c(dim.fixed, dim.free))) {  #There's a problem
    if(any(dim.fixed)) {
      msg = c(msg, paste("fixed", par.names[dim.fixed], "dims are incorrect. Dim 3 should be 1 or the number of time steps in the data.\n"))
    }
    if(any(dim.free)) {
      msg = c(msg, paste("free", par.names[dim.fixed], "dims are incorrect. Dim 3 should be 1 or the number of time steps in the data.\n"))
    }
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }

  ###########################
  # Check that free has column names since these are the parameter names
  ###########################
  no.colnames.free=unlist(lapply(lapply(free[par.names],colnames),is.null)) & unlist(lapply(free[par.names],function(x){dim(x)[2]}))!=0
  if(any(no.colnames.free)) {
    msg = c(msg, paste("free", par.names[no.colnames.free], "is missing column names.\n"))
  }

  ###########################
  # Check data and missing values consistency if data present
  # as.numeric(NA) is the missing value
  ###########################
  if(!is.numeric(MODELobj$data)) msg = paste(msg, "Data must be numeric. \n")
  for( bad.val in c(NA, NaN, Inf, -Inf)){
    if(!identical(bad.val, as.numeric(NA)) && ( bad.val %in% MODELobj$data ) ){  
      msg = c(msg, paste("Data cannot have ", bad.val,". \n",sep="")) }
  }

  ###########################
  # Y.names against rownames of data.  They should be identical; otherwise something got scrambled
  ###########################
  if(!identical(attr(MODELobj,"Y.names"),rownames(MODELobj$data)) ){  
    msg = c(msg, paste("The rownames of the data and the attribute Y.names don't match.\n",sep="")) }
  
  ###########################
  # Check tinitx; must be 0 or 1
  ###########################
  if( !(MODELobj$tinitx %in% c(0,1)) ){ msg = c(msg, "tinitx (t associated with initial x) must be 0 or 1.\n") }
  
  if(!is.null(msg)){  #next test won't work so stop now
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }

  ###########################
  # Last Check that fixed, free, par.names are complete and consistent
  # This is form dependent so the MARSS.form file needs to include a is.marssMODEL_form() function
  ###########################
  form=attr(MODELobj, "form")
  is.marssMODEL.fun = paste("is.marssMODEL_", form[1], sep="")
  tmp=try(exists(is.marssMODEL.fun, mode="function"),silent=TRUE)
  if(isTRUE(tmp)){
    #the is.marssMODEL_form function runs tests and then returns msgs
    msg=eval(call(is.marssMODEL.fun, MODELobj, method=method))
  }else{ 
    msg=c(msg, paste("No is.marssMODEL_", form[1], " is available to test the model object.\n", sep=""))
  }

  ###########################
  # Check diffuse; must be TRUE or FALSE
  ###########################
  if( !(MODELobj$diffuse %in% c(FALSE, TRUE)) ){ msg = c(msg, "diffuse must be TRUE or FALSE.\n") }
  
  if(length(msg) == 0){ return(TRUE)
  }else{
    msg=c("\nErrors were caught in is.marssMODEL()\n", msg)
    return(msg)
  }
   
}

