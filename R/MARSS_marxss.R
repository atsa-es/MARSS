###################################################################################
# marxss form definition file
# MARSS.marxss takes MARSS() call input 
# and returns the MARSS.call list with 2 model objects added
#    marssMODEL(form=marss) in the  $marss element
#    marssMODEL(form=marxss) in the $model element
# marss_to_marxss converts marssMODEL(form=marss) to marssMODEL(form=marxss)
# marxss_to_marss converts marssMODEL(form=marxss) to marssMODEL(form=marss)
# plus the following methods helper functions
# print_marxss, coef_marxss, predict_marxss, describe_marxss, MARSSinits_marxss
###################################################################################

###################################################################################
# Coerce model list input in MARSS() call to marxss model object (class marxss) put in $alt.forms
# Coerce marxss model into marss model and put in $marss
# using form = MARXSS
# x(t)=B(t)x(t-1) + U(t) + C(t)c(t) + w(t), W~MVN(0,Q)
# y(t)=Z(t)x(t) + A(t) + D(t)d(t) + v(t), V~MVN(0,R)
# x(t0) = x0 + l, L ~ MVN(0,V0)
# produces model object with fixed, free, tinitx, diffuse and data
# rhs is specified by fixed and free; lhs is specified with data (really should be list(x=, y=))
# attributes of model object has model.dims, X.names, form, equation

# The conversion functions has 3 parts
# Part 1 Set up the defaults and allowed structures
# Part 2 Do the conversion of model list to a marxss object
# Part 3 Do the conversion of marxss object to marss object

###################################################################################
MARSS.marxss=function(MARSS.call){
  #load needed package globals
  common.allowed.in.MARSS.call=get("common.allowed.in.MARSS.call", envir=pkg_globals)
  
  #Part 1 Set up defaults and check that what the user passed in is allowed
  
  #Check that no args were passed into MARSS that are not allowed
  marxss.allowed.in.MARSS.call=c("model")
  allowed.in.call=c(marxss.allowed.in.MARSS.call,common.allowed.in.MARSS.call)
  if(any(!(names(MARSS.call) %in% allowed.in.call))){
    bad.names=names(MARSS.call)[!(names(MARSS.call) %in% allowed.in.call)]
    msg=paste("Argument ", paste(bad.names, collapse=", "),"  not allowed MARSS call for form ", MARSS.call$form, ". See ?MARSS.marxss\n", sep="")
    cat("\n","Errors were caught in MARSS.marxss \n", msg, sep="")
    stop("Stopped in MARSS.marxss() due to problem(s) with model specification.\n", call.=FALSE)
  }
  
  # 1 Check for form dependent user inputs for method and reset defaults for inits and control if desired
  
  # 2 Specify the text shortcuts and whether factors or matrices can be passed in
  #   The names in the allowed list do not need to be A, B, Q .... as used in form=marss object
  #   Other names can be used if you want the user to use those names; then in the MARSS.form function
  #   you convert the user passed in names into the form marss names with the A, B, Q, R, ... names
  #   checkModelList() will check what the user passes in against these allowed values, so
  #   so you need to make sure each name in model.defaults has a model.allowed value here 
  model.allowed = list(
    A=c("scaling", "unconstrained", "unequal", "equal", "zero"),
    D=c("unconstrained", "equal", "unequal", "zero","diagonal and unequal","diagonal and equal","zero"),
    B=c("identity", "zero", "unconstrained", "unequal", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    Q=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    R=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    U=c("unconstrained", "equal", "unequal", "zero"),
    C=c("unconstrained", "equal", "unequal", "zero","diagonal and unequal","diagonal and equal","zero"),
    x0=c("unconstrained", "equal", "unequal", "zero","diagonal and unequal","diagonal and equal"),
    V0=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    Z=c("identity","unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "onestate"),
    G=c("identity","zero"),
    H=c("identity","zero"),
    L=c("identity","zero"),
    c=c("zero"),
    d=c("zero"),
    tinitx=c(0,1),
    diffuse=c(TRUE,FALSE),
    factors = c("Z"),
    matrices = c("A","B","Q","R","U","x0","Z","V0","D","C","d","c","G","H","L")
  )
  
  #model.defaults is form dependent so you must specify it
  model.defaults =list(Z="identity", A="scaling", R="diagonal and equal", B="identity", U="unconstrained", 
                       Q="diagonal and unequal", x0="unconstrained", V0="zero", D="zero", d=matrix(0,1,1), 
                       C="zero", c=matrix(0,1,1), G="identity", H="identity", L="identity", tinitx=0, diffuse=FALSE)
  
  if(!is.null(MARSS.call[["model"]][["c"]]))
    if(!identical(MARSS.call$model$c, "zero") & !all(MARSS.call$model$c==0)) model.defaults$C="unconstrained"
  if(!is.null(MARSS.call[["model"]][["d"]]))
    if(!identical(MARSS.call$model$d, "zero") & !all(MARSS.call$model$d==0)) model.defaults$D="unconstrained"
  
  #This checks that what user passed in model list can be interpreted and converted to form marss
  #if no errors, it updates the model list by filling in missing elements with the defaults
  MARSS.call$model=checkModelList( MARSS.call$model, model.defaults, model.allowed )
  
  # Part 2 Convert the model list to a marssMODEL object, form=marxss
  ## set up fixed and free elements
  fixed = free = list()
  
  model = MARSS.call[["model"]]
  model.elem = c("Z","A","R","B","U","Q","x0","V0","D","d","C","c","G","H","L") 
  dat = MARSS.call[["data"]]
  if(is.vector(dat)) dat=matrix(dat,1,length(dat))
  n = dim(dat)[1]; TT = dim(dat)[2]
  if(is.null(rownames(dat))){
    Y.names = paste("Y",seq(1, n),sep="") #paste(seq(1, n), sep="") 
    rownames(dat)=Y.names
  }else{ 
    Y.names = rownames(dat)
    if(any(duplicated(Y.names))){
      for(i in Y.names[duplicated(Y.names)]){
        loc = (Y.names==i)
        nn = sum(loc)
        Y.names[loc]=paste(i,"-",1:nn,sep="")
        rownames(dat)[loc] = Y.names[loc]
        MARSS.call[["data"]] = dat
      }
    }
  }
  
  ## Set m based on Z specification IF Z was specified; errors will be reported later if m conflicts with other parameters
  m = NA
  if (identical(model$Z, "unconstrained")) m = n
  if (identical(model$Z, "equalvarcov")) m = n
  if (identical(model$Z, "diagonal and equal")) m = n
  if (identical(model$Z, "diagonal and unequal")) m = n
  if (identical(model$Z, "onestate")) m = 1
  if (identical(model$Z, "identity")) m = n
  if (is.factor(model$Z)) m = length(levels(model$Z)) 
  if (is.array(model$Z)) m = dim(model$Z)[2] 
  
  X.names=NULL
  if(!(is.null(model[["X.names"]]))) X.names = model[["X.names"]]
    if(is.null(X.names) & identical(model$Z,"identity"))
      X.names=paste("X.",Y.names,sep="")
  if( is.null(X.names) && is.array(model$Z)){
    if(length(dim(model$Z))==3)
      if(dim(model$Z)[3]==1)
        if(is.design(model$Z) & !is.null(colnames(model$Z)))
          X.names=colnames(model$Z)
  }
  if( is.null(X.names) && is.factor(model$Z)){
    X.names=unique(as.character(model$Z))
  }
  if(is.null(X.names) && is.matrix(model$Z))
    if(is.design(model$Z) && !is.null(colnames(model$Z))) X.names=colnames(model$Z)
  if(is.null(X.names) && is.matrix(model$Z) && is.identity(model$Z)) 
      X.names = paste("X.",Y.names,sep="")
  if(is.null(X.names)) X.names = paste("X",seq(1, m),sep="") #paste(seq(1, m),sep="")  #
  
  ## Set c1 based on C specification if C specified with a particular shape
  ## error checking later will complain if C and c (or D and d) conflict
  if (is.array(model$C)) c1 = dim(model$C)[2] else c1 = 1
  if (is.array(model$D)) d1 = dim(model$D)[2] else d1 = 1

  ## Set based on G and H specification
  ## error checking later will complain if conflict
  if (is.array(model$G)) g1 = dim(model$G)[2] else g1 = m
  if (is.array(model$H)) h1 = dim(model$H)[2] else h1 = n
  if (is.array(model$L)) l1 = dim(model$L)[2] else l1 = m
  
  for(el in c("c","d")){
    thedim=get(paste(el,"1",sep=""))
    if(identical(model[[el]], "zero")){ model[[el]]=matrix(0,thedim,1) }
    #if(!any(is.na(model[[el]])) & all(model[[el]]==0)) model[[toupper(el)]]="zero"
    if(is.vector(model[[el]])) model[[el]]=matrix(model[[el]],1,length(model[[el]]))
  }
  #Now set c1 and d1 based on c and d, which should now be a matrix of some sort.  This ensures that c1 and d1 are set
  #if c and C or d and D conflict, this will be picked up the error checking because dims of C or D and model.dims won't match
  c1=dim(model$c)[1]; d1=dim(model$d)[1]
  #3rd dim of params are set to 1 and will be reset to correct value at end
  model.dims = list(data=c(n,TT),x=c(m,TT),y=c(n,TT),w=c(m,TT),v=c(n,TT),Z=c(n,m,1),U=c(m,1,1),
                    A=c(n,1,1),B=c(m,m,1),Q=c(g1,g1,1),R=c(h1,h1,1),x0=c(m,1,1),V0=c(l1,l1,1),
                    D=c(n,d1,1), d=c(d1,1,1), C=c(m,c1,1), c=c(c1,1,1), G=c(m,g1,1), H=c(n,h1,1), L=c(m,l1,1))
  
  ## Error-checking section that is specific to marxss form
  # Note most error checking happens in checkMARSSInputs, checkModelList, and is.marssMLE
  #If a model elements passed in as a factors, make sure it is the correct length otherwise construction of marssMODEL object will break
  problem = FALSE
  msg=NULL
  ## Check model structures that are passed in as factors
  correct.factor.len = list(Z=n,A=n,R=n,B=m,U=m,Q=m,x0=m,V0=m)  
  for (el in model.elem) {
    #if a factor then it needs to have the correct length otherwise construction of marssMODEL object will break and no NAs
    if( is.factor(model[[el]]) ) {
      if(length(model[[el]]) != correct.factor.len[[el]]) {
        problem=TRUE
        msg = c(msg, paste(" The model$", el, " is being passed as a factor and should be length ", correct.factor.len[[el]], " based on data dims. It's not. See help file.\n", sep="")) }
      if(NA.in.fac <- NA %in% model[[el]]) {
        problem=TRUE
        msg = c(msg, paste(" NAs are not allowed in model factor for ", el, ". See help file.\n", sep="")) }
    }  #is factor  
  } # end for (el in model.elem)
  
  #if el == Z then factor needs to have m levels
  if( is.factor(model$Z) ) {
    if(length(levels(model$Z)) != m) {
      problem=TRUE
      msg=c(msg," When Z is a factor, the number of levels must equal the number of state processes (m).\n")
    } }
  
  #Check that if A is scaling, then Z spec must lead to a design matrix
  if(identical(model$A,"scaling")){
    if(is.array(model$Z) & length(dim(model$Z))==3)
      if(dim(model$Z)[3]!=1 && !is.design(model$Z, zero.cols.ok = TRUE)) { #if it is a matrix
        problem=TRUE
        msg = c(msg, " If A is scaling(the default), then Z must be a time-constant design matrix:(0,1) and rowsums=1.\nYou can construct a scaling A matrix and pass that in.\n")
      }
    if((!is.array(model$Z) & !is.factor(model$Z))) #if it is a string
      if(!(model$Z %in% c("onestate","identity"))){
        problem=TRUE
        msg = c(msg, " If A is scaling(the default), then Z must be a time-constant design matrix:(0,1) and rowsums=1.\n")
      }
    if(is.matrix(model$Z) && !is.design(model$Z, zero.cols.ok = TRUE)) { #if it is a matrix, won't be array due to first test
      problem=TRUE
      msg = c(msg, " If A is scaling(the default), then Z must be a time-constant design matrix:(0,1) and rowsums=1.\n")
    }
  }  
  if(is.array(model$x0) & length(dim(model$x0))==3)
    if(dim(model$x0)[3]!=1){
      problem=TRUE
      msg = c(msg, " x0 cannot be time-varying. If x0 in model arg is 3D, the 3rd dim must equal 1.\n")
    } 
  if(is.array(model$V0) & length(dim(model$V0))==3)
    if(dim(model$V0)[3]!=1){
      problem=TRUE
      msg = c(msg, " V0 cannot be time-varying. If V0 in model arg is 3D, the 3rd dim must equal 1.\n")
    } 
  if(is.array(model$L) & length(dim(model$L))==3)
    if(dim(model$L)[3]!=1){
      problem=TRUE
      msg = c(msg, " V0 and thus L cannot be time-varying. If L in model arg is 3D, the 3rd dim must equal 1.\n")
    }
  #if C is diagonal and equal or diagonal and unequal, then d1=m
  if(identical(model$C,"diagonal and equal") || identical(model$C,"diagonal and unequal"))
    if(c1!=m){
      problem=TRUE
      msg = c(msg, " If C is diagonal, it must be square and c must be m x 1.\n")
    }
  #if D is diagonal and equal or diagonal and unequal, then d1=n
  if(identical(model$D,"diagonal and equal") || identical(model$D,"diagonal and unequal"))
    if(d1!=n){
      problem=TRUE
      msg = c(msg, " If D is diagonal, it must be square and d must be n x 1.\n")
    }
  #if c and d can't have any NAs or Infs
  for(el in c("c","d"))
    if(any(is.na(model[[el]])) || !is.numeric(model[[el]]) || any(is.infinite(model[[el]]))){
      problem=TRUE
      msg = c(msg, paste(" ",el,"must be numeric and have no NAs, NaNs, or Infs.\n"))
    }
  #c and d must be a 2D matrix and 2nd dim must be 1 or TT or a 3D matrix with 2nd dim = 1
  for(el in c("c","d")){
    if( length(dim(model[[el]]))!=2 & length(dim(model[[el]]))!=3){
      problem=TRUE
      msg = c(msg, paste(" ",el,"must be a 2D matrix with time in 2nd dim or 3D with time in 3rd dim.\n"))
    }else{
      if(length(dim(model[[el]]))==3){ #3D
        if( dim(model[[el]])[2]!=1 ){
          problem=TRUE
          msg = c(msg, paste(" if",el,"is 3D, 2nd dim must be 1 and time is in 3rd dim.\n"))
        }
        if( !(dim(model[[el]])[3] == 1 || dim(model[[el]])[3] == TT) ){
          problem=TRUE
          msg = c(msg, paste(" if",el,"is 3D, 3rd dim equal to 1 or T (length of data).\n"))
        }
      }else{ #is matrix
        if( !(dim(model[[el]])[2] == 1 || dim(model[[el]])[2] == TT) ){
        problem=TRUE
        msg = c(msg, paste(" ",el,"must be a 2D matrix with 2nd dim equal to 1 or T (length of data).\n"))
        }
      }
    }
  }
  
  #If there are problems
  if(problem)  {
    cat("\n","Errors were caught in MARSS.marxss \n", msg, sep="")
    stop("Stopped in MARSS.marxss() due to problem(s) with model specification.\n", call.=FALSE)
  }
  #end of error section  
  
  #If 2D matrix, change c and d to array so that it can be handled by the normal fixed/free constructions
  for(el in c("c","d")){
    if(is.matrix(model[[el]])){
    row.names=rownames(model[[el]])
    model[[el]]=array(model[[el]],dim=c(dim(model[[el]])[1],1,dim(model[[el]])[2]))
    rownames(model[[el]])=row.names
    }
  }
  
  
  ##Translate the text shortcuts into a marssMODEL object 
  ## Translate the model structure names (shortcuts) into fixed and free
  ## fixed is a dim(1)*dim(2) X 1 vector of the fixed (intercepts) values
  ## free is a dim(1)*dim(2) X p vector of the free (betas) values for the p estimated elements
  
  model.elem = c("Z","A","R","B","U","Q","x0","V0","D","C","d","c","G","H","L")
  if(which(model.elem=="Z")>which(model.elem=="A")) model.elem=rev(model.elem) #Z must go first
  
  tmp=list()
  for(el in model.elem) {
    if(MARSS.call$silent==2) cat(paste("  Building fixed and free matrices for ",el,".\n",sep=""))
    tmp[[el]]="not assigned"
    if(el=="Z" & is.factor(model$Z)) { 
      tmp[[el]] = matrix(0,model.dims$Z[1], model.dims$Z[2])  
      for(i in X.names) tmp[[el]][which(model$Z==i), which(as.vector(X.names)==i)] = 1
    }
    if(el=="Z" & identical(model$Z,"onestate")) {   #m=1
      tmp[[el]] = matrix(1, n, 0)
    }
    if( identical(model[[el]],"identity") ) { 
      tmp[[el]] = diag(1,model.dims[[el]][1])
    }
    if(identical(model[[el]],"diagonal and equal")) {
      tmp[[el]] = array(list(0),dim=c(model.dims[[el]][1],model.dims[[el]][2]))
      diag(tmp[[el]])="diag" #paste(el,"(diag)",sep="")
      if(length(tmp[[el]])==1) tmp[[el]][1,1]=el
    }
    if(identical(model[[el]],"diagonal and unequal")) {
      tmp[[el]] = array(list(0),dim=c(model.dims[[el]][1],model.dims[[el]][2]))
      dim.mat = model.dims[[el]][1]
      el.labs=as.character(1:dim.mat)
      if(el %in% c("V0","Q","B")) el.labs=X.names
      if(el %in% c("Z","R")) el.labs=Y.names
      diag(tmp[[el]])=paste("(",el.labs,",", el.labs,")",sep="") #paste(el,"(",as.character(1:dim.mat),",",as.character(1:dim.mat),")",sep="")
      if(length(tmp[[el]])==1) tmp[[el]][1,1]=el
    }
    if(identical(model[[el]],"unconstrained") || identical(model[[el]],"unequal")){
      tmp[[el]]=array(NA,dim=c(model.dims[[el]][1],model.dims[[el]][2])) 
      if(el %in% c("Q","R","V0")){  #variance-covariance matrices
        dim.mat = model.dims[[el]][1]
#         for(i in 1:dim.mat){
#           for(j in 1:dim.mat) tmp[[el]][i,j]=tmp[[el]][j,i]=paste("(",i,",",j,")",sep="") #paste(el,"(",i,",",j,")",sep="")
#         }
          tmp[[el]]=matrix(paste("(",rep(1:dim.mat,dim.mat),",",rep(1:dim.mat,each=dim.mat),")",sep=""),dim.mat,dim.mat)
          tmp[[el]][upper.tri(tmp[[el]])] = t(tmp[[el]])[upper.tri(t(tmp[[el]]))]  
      }else{ #not var-cov matrix
        row.name=1:model.dims[[el]][1]
        col.name=1:model.dims[[el]][2]
        if(el %in% c("U","x0")) row.name=X.names
        if(el == "A") row.name=Y.names
        if(el %in% c("C","D")){
          if(el=="C") row.name=X.names
          if(el=="D") row.name=Y.names
          if(!is.null(rownames(model[[tolower(el)]]))) col.name=rownames(model[[tolower(el)]])
        }
#         for(i in 1:model.dims[[el]][1]){
#           for(j in 1:model.dims[[el]][2]){
#             if(model.dims[[el]][2]>1) tmp[[el]][i,j]=paste("(",row.name[i],",",col.name[j],")",sep="") #paste(el,"(",row.name[i],",",col.name[j],")",sep="")
#             else tmp[[el]][i,j]=paste(row.name[i],sep=",") #paste(el,row.name[i],sep=",")
#           }
#         }
        if(model.dims[[el]][2]>1){
          tmp[[el]]=matrix(paste("(",rep(row.name,model.dims[[el]][2]),",",rep(col.name,each=model.dims[[el]][1]),")",sep=""),model.dims[[el]][1],model.dims[[el]][2])
        }else{
          tmp[[el]]=matrix(row.name,model.dims[[el]][1],1)
        }
      }
      if(length(tmp[[el]])==1) tmp[[el]][1,1]=el 
    } #unconstrained
    if(identical(model[[el]],"equalvarcov")) {
      tmp[[el]]=array("offdiag",dim=c(model.dims[[el]][1],model.dims[[el]][2])) #array(paste(el,"(offdiag)",sep=""),dim=model.dims[[el]])
      diag(tmp[[el]])="diag" #paste(el,"(diag)",sep="")
      if(length(tmp[[el]])==1) tmp[[el]][1,1]=el
    }
    if(identical(model[[el]],"equal")) { 
      tmp[[el]]=array("1",dim=c(model.dims[[el]][1],model.dims[[el]][2])) #array(el,dim=model.dims[[el]])
    }
    if(identical(model[[el]],"zero")) { 
      tmp[[el]]=array(0,dim=c(model.dims[[el]][1],model.dims[[el]][2]))
    }
    if(is.array(model[[el]])) {
      tmp[[el]]=model[[el]]
    }
    if(el=="A" & identical(model[[el]], "scaling")) {  #check above ensures that Z is design and time-invariant
      ## Construct A from fixed Z matrix
      tmp[[el]] = matrix(list(),model.dims$A[1],model.dims$A[2])
      tmp[[el]][,1]=Y.names  
      for(i in 1:m) {
        nonzeroZ=tmp$Z[,i]!=0
        if(any(nonzeroZ)) tmp[[el]][min(which(nonzeroZ)), 1] = 0
      }
    }     
    if(identical(tmp[[el]],"not assigned")) stop(paste("Stopped in MARSS.marxss(): tmp was not assigned for ",el,".\n",sep=""))
    tmpconst=convert.model.mat(tmp[[el]])
    free[[el]] = tmpconst$free
    fixed[[el]] = tmpconst$fixed 
    
    #set the last dim of the model.dims since it was at a temp value to start
    model.dims[[el]][3]=max(dim(free[[el]])[3],dim(fixed[[el]])[3])
  }
  
  #save the row names for the inputs by setting in fixed
  for(el in c("c","d")){
    if(is.null(rownames(tmp[[el]]))) rownames(tmp[[el]])=paste(el,seq(1,dim(tmp[[el]])[1]),sep="")
    rownames(fixed[[el]])=rownames(tmp[[el]])
  }
  
  
  #Set the marssMODEL form marxss
  #This is the f+Dp form for the MARXSS model used for user displays, printing and such
  marxss_object = list(fixed=fixed, free=free, data=dat, tinitx=model$tinitx, diffuse=model$diffuse)
  #set the attributes
  class(marxss_object) <- "marssMODEL"
  attr(marxss_object, "obj.elements") <- c("fixed","free","data","tinitx","diffuse")
  attr(marxss_object, "form") <- "marxss"
  attr(marxss_object, "model.dims") <- model.dims
  #par.names are what needs to be in fixed/free pair
  attr(marxss_object, "par.names") <- c("Z","A","R","B","U","Q","x0","V0","D","C","d","c","G","H","L")
  attr(marxss_object, "X.names") <- X.names
  attr(marxss_object, "Y.names") <- Y.names
  attr(marxss_object, "equation") <- "x_{t}=B_{t}*x_{t-1}+U_{t}+C_{t}*c_{t}+G{t}*w_{t}; w_{t}~MVN(0,Q_{t})\ny_{t}=Z_{t}*x_{t}+A_{t}+D_{t}*d_{t}+H{t}*v_{t}; v_{t}~MVN(0,R_{t})"

  #Change alldefaults global to match the form
  #first load the defaults
  alldefaults <- get("alldefaults", envir=pkg_globals)
  alldefaults[[MARSS.call$method]][["inits"]][["C"]] <- 0
  alldefaults[[MARSS.call$method]][["inits"]][["D"]] <- 0
  #c and d and G and H inits won't be used but assigned defaults so users can pass in inits as coef(fit)
  alldefaults[[MARSS.call$method]][["inits"]][["c"]] <- 0
  alldefaults[[MARSS.call$method]][["inits"]][["d"]] <- 0
  alldefaults[[MARSS.call$method]][["inits"]][["G"]] <- 0
  alldefaults[[MARSS.call$method]][["inits"]][["H"]] <- 0
  alldefaults[[MARSS.call$method]][["inits"]][["L"]] <- 0
  assign("alldefaults", alldefaults, pkg_globals)
  
  ## Check that the marssMODEL object output by MARSS.form() is ok since marxss_to_marss will go south otherwise
  if(MARSS.call$silent==2) cat(paste("  Running is.marssMODEL on the marxss model.\n",sep=""))
  tmp <- is.marssMODEL(marxss_object, method=MARSS.call$method)
  if(!isTRUE(tmp)) {
    cat(tmp) 
    stop("Stopped in MARSS.marxss() due to problem(s) with model specification.", call.=FALSE)
  }
  
  #Put the marxss model into $model since model holds the model in the 'form' form
  MARSS.call$model <- marxss_object
  
  ## Create marssMODEL(form=marss) object added to call
  #when called with a marssMODEL object (as here), marxss_to_marss returns a marssMODEL object
  MARSS.call$marss <- marxss_to_marss(marxss_object)
  
  ## Return MARSS call list with $marss and $model added
  MARSS.call
}

marss_to_marxss=function(x, C.and.D.are.zero=FALSE){
  if(!(inherits(x, "marssMLE") || inherits(x, "marssMODEL"))) stop("Stopped in marss_to_marxss(): this function needs a marssMODEL or marssMLE object")
  #This function returns a MLE object where the model and par parts of the MLE object are in marxss form for printing purposes.
  #This function needs a marxss marssMODEL object and will break otherwise
  #You cannot back construct a marxss from the marss model
  #The function will also work is x is a model object, but then it just returns marxss.marssMODEL
  #written this way so it doesn't crash if x happens to be a marssMODEL in case I later dynamically write function names/calls
  
  if(inherits(x, "marssMODEL")){
    marss.model=x
    if(!("marss" %in% attr(marss.model,"form"))) stop("Stopped in marss_to_marxss(): this function requires a marssMODEL object in marss form.\n",call.=FALSE)
  }else{ 
    marss.model=x[["marss"]]
    if(!("marss" %in% attr(marss.model,"form"))) stop("Stopped in marss_to_marxss(): this function requires a marssMLE object with element $marss a marssMODEL in marss form.\n",call.=FALSE)
  }
  
  if(!C.and.D.are.zero){
    if(inherits(x, "marssMODEL")){ 
      stop("Stopped in marss_to_marxss(: function was called with a marss model object instead of MLE object, so needs a marxss model passed in.")
    }else{ 
      marxss.model=x[["model"]] #should be model since we want the marxss form
      if(any(is.null(attr(marxss.model,"form")),!("marxss" %in% attr(marxss.model,"form")))){
        stop("Stopped in marss_to_marxss(: function was called with a MLE object, so needs the $model element to be form marxss.")
      }
    }
  }else{ #C and D are zero so we can construct a marxss object
    marxss.model = marss.model #use the marss model as a template
    marxss.dims=attr(marss.model,"model.dims")
    marxss.model[["fixed"]][["C"]]=array(0,dim=c(marxss.dims[["x"]][1],1,1))
    marxss.model[["fixed"]][["D"]]=array(0,dim=c(marxss.dims[["y"]][1],1,1))
    marxss.model[["free"]][["C"]]=array(0,dim=c(marxss.dims[["x"]][1],0,1))
    marxss.model[["free"]][["D"]]=array(0,dim=c(marxss.dims[["y"]][1],0,1))
    marxss.model[["fixed"]][["c"]]=array(0,dim=c(1,1,1))
    marxss.model[["fixed"]][["d"]]=array(0,dim=c(1,1,1))
    marxss.model[["free"]][["c"]]=array(0,dim=c(1,0,1))
    marxss.model[["free"]][["d"]]=array(0,dim=c(1,0,1))
    marxss.dims[["C"]]=c(marxss.dims[["x"]][1],1,1)
    marxss.dims[["D"]]=c(marxss.dims[["y"]][1],1,1)
    marxss.dims[["c"]]=c(1,1,1)
    marxss.dims[["d"]]=c(1,1,1)
    #reset the attributes to marxss form
    #obj.elements, X.names and Y.names stay the same
    attr(marxss.model, "form")="marxss"
    attr(marxss.model, "model.dims")=marxss.dims
    #par.names are what needs to be in fixed/free pair; order is important
    attr(marxss.model, "par.names")=c("Z","A","R","B","U","Q","x0","V0","D","C","d","c","G","H","L")
    attr(marxss.model, "equation")="x_{t}=B_{t}*x_{t-1}+U_{t}+C_{t}*c_{t}+G{t}*w_{t}; w_{t}~MVN(0,Q_{t})\ny_{t}=Z_{t}*x_{t}+A_{t}+D_{t}*d_{t}+H{t}*v_{t}; v_{t}~MVN(0,R_{t})"
  }
  if(inherits(x, "marssMODEL")) return(marxss.model) #in marxss form

  x[["model"]]=marxss.model #now in marxss form
  marxss.dims = attr(marxss.model, "model.dims")
  
  #got here, so class of x is marssMLE
    for(val in c("par","start","par.se","par.bias","par.upCI","par.lowCI")){
      if(!is.null(x[[val]])){
          tmp.dim=dim(marxss.model$free$C)[2] #how many C parameters
        if(tmp.dim==0){
          x[[val]][["C"]] = matrix(0,0,1)
        }else{
          #because marss.U is [marxss.C marxss.U]
          x[[val]][["C"]] = x[[val]][["U"]][1:tmp.dim,, drop=FALSE]
          x[[val]][["U"]] = x[[val]][["U"]][-(1:tmp.dim),, drop=FALSE]
        }
          tmp.dim=dim(marxss.model$free$D)[2] #how many D parameters
        if(tmp.dim==0){
          x[[val]][["D"]] = matrix(0,0,1)
        }else{
          #because marss.A is [marxss.D marxss.A]
          x[[val]][["D"]] = x[[val]][["A"]][1:tmp.dim,, drop=FALSE]
          x[[val]][["A"]] = x[[val]][["A"]][-(1:tmp.dim),, drop=FALSE]
        }
        for(el in c("c","d")){
          x[[val]][[el]] = matrix(0,0,1) #because c and d are inputs not estimated
        }
      } #not is null val
    } #for val in par, start

    #returning a marssMLE object where the par, start etc are in marxss form
  # $model is in marxss and $marss stays the same
  return(x) 
}

marxss_to_marss=function(x, only.par=FALSE){
  #x is a marssMODEL of form marxss
  #this will create a marss model object (if !only.par) and a par list in form marss
  #if only.par=TRUE then only the par element is changed and marss is used for the marss object
  
  #hold on to this since x will be changing and need to know what to return
  class.x=class(x)
  if(!(class.x %in% c("marssMODEL","marssMLE"))){
    stop("Stopped in marss_to_marxss(): x$model must be a marssMODEL or marssMLE.",call.=FALSE)
  }
  
  #check form if user passed in marssMODEL
  if(class.x=="marssMODEL"){
    if(!("marxss" %in% attr(x, "form"))) stop("Stopped in marss_to_marxss(): this function requires a marssMODEL object in marxss form.")
  }
  
  if(class.x=="marssMLE"){ #Then set the par elements to correspond to marss if they are in marxss form
    #check that the model element they passed in is marxss
    if(!("marxss" %in% attr(x[["model"]], "form"))) stop("Stopped in marss_to_marxss(): x$model must be in marxss form.",call.=FALSE)
    
    x.marss = list()
    
    #this will go through x and reset any par-like obj in marxss form
    for(val in c("par","start","par.se","par.bias","par.upCI","par.lowCI")){
      if(!is.null(x[[val]])){
          #because marss.U is [marxss.C marxss.U]
          x.marss[[val]][["U"]] = rbind(x[[val]][["C"]],x[[val]][["U"]])
          #because marss.A is [marxss.D marxss.A]
          x.marss[[val]][["A"]] = rbind(x[[val]][["D"]],x[[val]][["A"]])
          #other elements are the same as for marxss
          for(el in c("R","Q","B","Z","x0","V0","G","H","L"))
            x.marss[[val]][[el]]=x[[val]][[el]]
          #reset x[[val]] so it only includes the marss elements
          x[[val]]=x.marss[[val]] #replace the x[[val]] with marss version
      } #not is null val
    } #for val in par, start
    
    marxss.model=x[["model"]]
  }else{ #end if x is marssMLE
    marxss.model = x
  }
  #marxss.model is a marssMODEL in marxss form and x is the original marssMLE object
  #if only.par updating was requested, return x
  if(class.x == "marssMLE" & only.par & !is.null(x[["marss"]])) return(x)
  
  #else construct a marss model object from marxss.model and put in $marss
  
  marxss.dims=attr(marxss.model,"model.dims")  
  fixed=marxss.model[["fixed"]]  #fixed and free will be modified, so these are holders not shortcuts
  free=marxss.model[["free"]]
  #marss.dims will be modified, so this is a holder not a shortcut
  marss.dims=marxss.dims
  n=marss.dims[["y"]][1]; m=marss.dims[["x"]][1]; TT=marss.dims[["x"]][2]
  #This step converts U+Cc into equivalent Uu and A+Dd into Aa
  #So U --> [C U] and u --> [c \\ 1]; A --> [D A] and a --> [d \\ 1]
  #This code adds fixed$u and fixed$a, and changes fixed$U and fixed$A; otherwise all stays the same
  for(el in c("C","D")){
    if(el=="C") el2="U" else el2="A"
    #if el all zero (fixed and all zero), it doesn't appear in the equation and U-->U and u-->1
    if(!is.fixed(marxss.model[["free"]][[el]]) | !all(sapply(marxss.model[["fixed"]][[el]],function(x){isTRUE(all.equal(x,0))}))){      
      #create fixed$u or fixed$a by add 1 to bottom of fixed$c or fixed$d
      #since fixed$c is a 3D array, we need to do this in an odd way:
      fixed[[tolower(el2)]]=array(apply(fixed[[tolower(el)]],3,rbind,1),dim=dim(fixed[[tolower(el)]])+c(1,0,0))
      #this fixed$u is just an input so we set the free to not estimated
      free[[tolower(el2)]]=array(0,dim=c(1,0,1)) #not estimated so 0 columns
      
      #next change U to [C U] and A to [D A]; need to figure out the dims since
      #C or U might be time-varying
      Tmax.fixed=max(dim(fixed[[el]])[3],dim(fixed[[el2]])[3])
      Tmax.free=max(dim(free[[el]])[3],dim(free[[el2]])[3])
      #dim. is c(dim 1 of C, dim 2 of C, dim 1 of U, dim 2 of U)
      dim.fixed=c(dim(fixed[[el]])[1],dim(fixed[[el]])[2],dim(fixed[[el2]])[1],dim(fixed[[el2]])[2])
      dim.free=c(dim(free[[el]])[1],dim(free[[el]])[2],dim(free[[el2]])[1],dim(free[[el2]])[2])
      #now that the dimensions are known, create and array holder for fixed$U and fixed$A
      tmp.fixed=array(NA,dim=c(dim.fixed[1]+dim.fixed[3],1,Tmax.fixed))
      tmp.free=array(0,dim=c(dim.free[1]+dim.free[3],dim.free[2]+dim.free[4],Tmax.free))
      #the first rows of fixed$U are fixed$C
      tmp.fixed[1:dim.fixed[1],,]=fixed[[el]]
      #the next rows are marxss.model$fixed$U
      tmp.fixed[(dim.fixed[1]+1):dim(tmp.fixed)[1],,]=fixed[[el2]]
      #Now create the new free$U.  If C estimated, free$C appears in the upper left
      if(dim.free[2]>0) tmp.free[1:dim.free[1],1:dim.free[2],]=free[[el]]
      #If U estimated, free$U appears in the lower right
      if(dim.free[4]>0) tmp.free[(dim.free[1]+1):dim(tmp.free)[1],(dim.free[2]+1):dim(tmp.free)[2],]=free[[el2]]
      #retain the column names (estimated parameter names)
      colnames(tmp.free)=c(colnames(free[[el]]),colnames(free[[el2]]))
      #assign the new fixed$U and free$U to fixed and free
      fixed[[el2]]=tmp.fixed
      free[[el2]]=tmp.free
      #settign U and A dims
      marss.dims[[el2]][2]=marxss.dims[[el]][2]+marxss.dims[[el2]][2]
    }else{  #Both C and U are all zero (fixed and all zero) 
      #so u is just 1 (a 1x1 matrix)
      fixed[[tolower(el2)]]=array(1,dim=c(1,1,1))
      free[[tolower(el2)]]=array(0,dim=c(1,0,1)) #not estimated
    }
  }
  
  #Now the fixed/free specify x=Bx+U(t)u(t)+w(t) and y=Zx+A(t)a(t)+v(t)
  #this part converts U(t)u(t) to U(t) and A(t)a(t) to A(t)
  #This requires making a vec(U(t)) that is specified by (rhs at end):
  #vec(Uu)=(t.u kron I_m)vec(U)=fixed+free*p=(t.u kron I_m)f+(t.u kron I_m)Dp
  #where f and D are fixed$U and free$U for U(t)u(t) form
  #the small case are inputs and the large case are estimated parameters  
  for(el in c("U","A")){
    #if c or a passed in. If not will be array(1,dim=c(1,1,1))
    #this if statement is just avoiding unneccesary code.  The math should still hold whether or not c is 1
    if(!identical(unname(fixed[[tolower(el)]]), array(1,dim=c(1,1,1)))){
      #hold onto fixed$U and free$U (not marxss.model$fixed and free but the new ones)
      free.orig=free[[el]]; fixed.orig=fixed[[el]]
      dim.free2=dim(free.orig)[2]; dim.free3=dim(free.orig)[3];  dim.fixed3=dim(fixed.orig)[3]; 
      dim.u.3=dim(fixed[[tolower(el)]])[3]
      Tmax=max(dim.free3, dim.fixed3, dim.u.3)
      #need the new marss U and A dims here which were defined above
      free[[el]]=array(0,dim=c(marss.dims[[el]][1],dim.free2,Tmax))
      colnames(free[[el]])=colnames(free.orig)
      fixed[[el]]=array(0,dim=c(marss.dims[[el]][1],1,Tmax))
      for(t in 1:Tmax){
        #the f and D of U (or A)
        f.t=sub3D(fixed.orig,t=min(t,dim.fixed3))
        d.t=sub3D(free.orig,t=min(t,dim.free3))
        #column vector of the u (or a) at time t
        ua.t=sub3D(fixed[[tolower(el)]],t=min(t,dim.u.3))  
        #vec(Uu)=(t.u kron I_m)vec(U)=(t.u kron I_m)f+(t.u kron I_m)Dp
        #again we want the new marss.dims defined above
        free[[el]][,,t]=(t(ua.t)%x%diag(1,marss.dims[[el]][1]))%*%d.t
        fixed[[el]][,,t]=(t(ua.t)%x%diag(1,marss.dims[[el]][1]))%*%f.t
      }
    }
  }
  marss.elem = c("Z","A","R","B","U","Q","x0","V0","G","H","L")
  free=free[marss.elem]
  fixed=fixed[marss.elem]
  dim3s=apply(rbind(unlist(lapply(free[marss.elem],function(x){dim(x)[3]})), unlist(lapply(fixed[marss.elem],function(x){dim(x)[3]}))),2,max)
  marss.dims = marxss.dims[!(names(marxss.dims) %in% c("C","c","D","d"))]
  marss.dims$U=c(m,1,dim3s[["U"]])
  marss.dims$A=c(n,1,dim3s[["A"]])
#   marss.dims = list(data=c(n,TT),x=c(m,TT),y=c(n,TT),w=c(m,TT),v=c(n,TT),
#                     Z=c(n,m,dim3s[["Z"]]),U=c(m,1,dim3s[["U"]]),A=c(n,1,dim3s[["A"]]),
#                     B=c(m,m,dim3s[["B"]]),Q=c(m,m,dim3s[["Q"]]),R=c(n,n,dim3s[["R"]]),
#                     x0=c(m,1,1),V0=c(m,m,1))
  
  ## Create the marss marssMODEL object
  marss.model = list(fixed=fixed, free=free, data=marxss.model[["data"]], tinitx=marxss.model[["tinitx"]], diffuse=marxss.model[["diffuse"]])
  #set the attributes that change
  class(marss.model) = "marssMODEL"
  attr(marss.model, "form")="marss"
  attr(marss.model, "obj.elements")=c("fixed","free","data","tinitx","diffuse")
  attr(marss.model, "model.dims")=marss.dims
  attr(marss.model, "par.names")=marss.elem
  attr(marss.model, "X.names")=attr(marxss.model,"X.names")
  attr(marss.model, "Y.names")=attr(marxss.model,"Y.names")
  attr(marss.model, "equation")="x_{t}=B_{t}*x_{t-1}+U_{t}+G{t}*w_{t}; w_{t}~MVN(0,Q_{t})\ny_{t}=Z_{t}*x_{t}+A_{t}+H{t}*v_{t}; v_{t}~MVN(0,R_{t})"
  if(class.x=="marssMODEL"){ return(marss.model) #marssMODEL of form marss
  }else{ 
    #class.x=marssMLE, then adds the marss element to the marssMLE object
    x[["marss"]]=marss.model
    return(x) #returning a marssMLE object where the par, start etc are in marss form and marss element is set
  }
}

#the par element of a marssMLE object is in form=marss.  Convert to form=marxss for printing
print_marxss = function(x){ return(marss_to_marxss(x)) } 

#the par element of a marssMLE object is in form=marss.  Convert to form=marxss for printing
coef_marxss = function(x){ return(marss_to_marxss(x)) } #this uses $model for marxss object

MARSSinits_marxss = function(MLEobj, inits){
  alldefaults=get("alldefaults", envir=pkg_globals)
  if(is.null(MLEobj[["model"]])){
    stop("Stopped in MARSSinits_marxss(): this function needs a marssMODEL in marxss form in $model",call.=FALSE)
  }else{
    if(!inherits(MLEobj[["model"]],"marssMODEL")) stop("Stopped in MARSSinits_marxss(): this function needs a marssMODEL in marxss form in $model",call.=FALSE)
    if(!("marxss" %in% attr(MLEobj[["model"]],"form"))) stop("Stopped in MARSSinits_marxss(): this function needs a marssMODEL in marxss form in $model",call.=FALSE)    
  }
    
  #B, Z, R, Q, x0 and V0 stay the same
  #U and A change
  #this function will return a U and A element for inits
  if(is.null(inits)) inits=list()
  elems=c("U","A","C","D")
  for(elem in elems){
    tmp.dim=dim(MLEobj$model$free[[elem]])[2] #how many estimated pars in marxss vers
    if(!is.null(inits[[elem]]) & !(tmp.dim==0)){ #tmp.dim==0 means no estimated
      if(!(length(inits[[elem]]) %in% c(tmp.dim,1))) 
        stop(paste("Stopped in MARSSinits_marxss(): ", elem," inits must be either a scalar (dim=NULL) or a matrix with 1 col and rows equal to the num of est values in ",elem,".",sep=""), call.=FALSE )
      if(tmp.dim!=0) inits[[elem]] = matrix(inits[[elem]],tmp.dim,1) else inits[[elem]]=matrix(0,0,1)
    }else{
      inits[[elem]] = matrix(alldefaults[[MLEobj$method]][["inits"]][[elem]],tmp.dim,1) }
  }
  inits$U = rbind(inits$C,inits$U) #yes, C on top
  inits$A = rbind(inits$D,inits$A) #yes, D on top
  return(inits)
}

predict_marxss = function(x, newdata, n.ahead, t.start){
  #x is a marssMLE object
  
  #This takes the newdata argument from a predict call and interprets the inputs in the context of the form of x$model
  #It will return a marssMODEL (form=marss) object ready for use in prediction
  #n.ahead is 1 by default; t.start is TT+1 by default
  #this makes tinitx=0, and uses E(x)_(t.start-1) as if t.start=1, init x is specified by x0T; V0 is V0T
  #if t.start!=1, init x is specified by xtT[,t.start-1] and VtT[,,t.start-1]
  
  #First get a par list in marxss form
  marxss.par = marss_to_marxss(x)[["par"]] #we only need par changed since marxss is in $model
  marxss.dims = attr(x[["model"]], "model.dims") #the marxss model dims
  marxss.model = x[["model"]]
  if( !("marxss" %in% attr(marxss.model,"form")) ) stop("Stopped in predict_marxss(): x$model needs to be in marxss form.", call.=FALSE)
  TT=marxss.dims[["y"]][2]
  form="marxss"
  allow.in=c("data","c","d") #allowed in newdata
  
  if(!(class(newdata)[1] %in% c("list","data.frame")))
    stop("Stopped predict_marxss(): newdata must be a list or dataframe.",call.=FALSE)
  
  #next interpret newdata;
  
  #if user passes in a data.frame, make an effort to interpret that
  #the following code makes sure newdata is a list with data, c and d with same 3rd dim (n.ahead)
  if(is.data.frame(newdata)){ #try to construct the list
    newdata.dataframe=newdata
    newdata=list()
    names.dataframe=names(newdata.dataframe)
    if(any(duplicated(names.dataframe))) stop("Stopped in predict_marxss(): the dataframe should not have any duplicated names",call.=FALSE)
    
    #first construct the data matrix from the dataframe
    Y.names=attr(marxss.model,"Y.names")
    #find any column names in the dataframe that match the rownames in the data matrix
    Y.match=match(Y.names,names.dataframe)
    if(any(!is.na(Y.match))){
      if(dim(newdata)[1]!=n.ahead){
        stop("Stopped in predict_marxss(): you have passed data in with newdata as a dataframe.\nThe number of rows must equal n.ahead in this case.",call.=FALSE)
      }
      cat(paste("Alert from predict_marxss(): y (data) are present in the dataframe and prediction will be conditioned on these values.\n",collapse=""))
      newdata[["data"]]=newdata.dataframe[Y.names[!is.na(Y.match)]]
    }
    #those that are missing are replaced with NA
    if(any(is.na(Y.match)))
      newdata[["data"]][Y.names[is.na(Y.match)]]=as.numeric(NA)
    #make sure the columns are in the same order as Y.names
    newdata[["data"]]=newdata$data[Y.names]
    newdata[["data"]]=t(as.matrix(newdata[["data"]])) #time across columns
    
    if(!(dim(newdata)[1]==n.ahead || dim(newdata)[1]!=1)){
      stop("Stopped in predict_marxss(): The number of rows in the newdata dataframe must be 1 or n.ahead.",call.=FALSE)
    }
    
    #Create the c and d matrices
    for(el in c("c","d")){
      is.zero = all(marxss.model[["fixed"]][[toupper(el)]]==0) & all(marxss.par[[toupper(el)]]==0)
      if(is.zero){ #if C or D is all zero then c or d is not needed; make it all 0
        newdata[[el]]=matrix(0,marxss.dims[[el]][1],1)
        rownames(newdata[[el]])=rownames(marxss.model[["fixed"]][[el]])
      }else{ #need c or d
        el.names=rownames(marxss.model[["fixed"]][[el]])        
        #find any column names in the dataframe that match the rownames in the el
        el.match=match(el.names,names.dataframe)
        #subset those columns in the dataframe that match the el names
        el.in.dataframe=el.names[!is.na(el.match)]
        newdata[[el]]=newdata.dataframe[el.in.dataframe]
        #make it into a matrix with time going across the columns like in $model
        newdata[[el]]=t(as.matrix(newdata[[el]]))
        #deal with any missing c or d rows
        if(!all(el.names %in% names.dataframe)){
          bad.names = el.names[!(el.names %in% names.dataframe)]
          if((n.ahead+t.start-1)>TT){
            #require that user passes in all the c and d inputs                       
            stop(c("Stopped in predict_marxss(): some of the ",el," inputs are missing: ",paste(bad.names,collapse=", ")),call.=FALSE)
          }else{ #replace missing names with values from the model
            if(marxss.dims[[el]][3]==1) t.el=1 else t.el=t.start:(t.start+n.ahead-1)
            tmp.el=matrix(marxss.model[["fixed"]][[el]][bad.names,1,t.el],length(bad.names),t.start+n.ahead-1)
            rownames(tmp.el)=bad.names
            newdata[[el]]=rbind(newdata[[el]],tmp.el)
          }
        }      
        #makesure the columns are in the same order as el.names
        newdata[[el]]=newdata[[el]][el.names,,drop=FALSE]      
        if(any(is.na(newdata[[el]]))){
          stop(paste("Stopped in predict_marxss(): there cannot be any NAs in the ",el," matrix.\n",sep=""),call.=FALSE)
        }
      }
    } #for el
  } #newdata is now a list with elements data, c, and d
  
  #if user passed in a list
  if(is.list(newdata)){
    for(el in c("c","d")){
      is.zero = all(marxss.model[["fixed"]][[toupper(el)]]==0) & all(marxss.par[[toupper(el)]]==0)
      if(is.zero){ #if C or D is all zero then c or d is not needed; make it all 0
        newdata[[el]]=matrix(0,marxss.dims[[el]][1],1)
        rownames(newdata[[el]])=rownames(marxss.model[["fixed"]][[el]])
      }else{ #need c or d
        if(!(el %in% names(newdata))){
          if((n.ahead+t.start-1)>TT){
            #require that user passes in all the c and d inputs                       
            stop(c("Stopped in predict_marxss(): Model has ", toupper(el), " but ", el," is missing from newdata\n and cannot be inferred since prediction extends beyond original dataset."),call.=FALSE)
          }else{ #replace missing names with values from the model
            if(marxss.dims[[el]][3]==1) t.el=1 else t.el=t.start:(t.start+n.ahead-1)
            newdata[[el]]=matrix(marxss.model[["fixed"]][[el]][,,t.el],marxss.dims[[el]][1],t.start+n.ahead-1)
            rownames(newdata[[el]])=rownames(marxss.model[["fixed"]][[el]])
          }
        }
        if(!is.matrix(newdata[[el]]))
          stop(c("Stopped in predict_marxss(): ", el," in newdata must be a matrix with time across the columns and n.ahead columns."),call.=FALSE)
        names.list=rownames(newdata[[el]])
        el.names=rownames(marxss.model[["fixed"]][[el]])
        #find any row names in the matrix that match the rownames in the el
        el.match=match(el.names,names.list)
        #subset those columns in the dataframe that match the el names
        el.in.matrix=el.names[!is.na(el.match)]
        newdata[[el]]=newdata[[el]][el.in.matrix,,drop=FALSE]
        if(!all(el.names %in% names.list)){
          bad.names = el.names[!(el.names %in% names.list)]
          if((n.ahead+t.start-1)>TT){
            #require that user passes in all the c and d inputs                       
            stop(c("Stopped in predict_marxss(): some of the ",el," inputs are missing: ",paste(bad.names,collapse=", ")," \nand cannot be inferred since prediction is beyond the end of the original dataset."),call.=FALSE)
          }else{ #replace missing names with values from the model
            if(marxss.dims[[el]][3]==1) t.el=1 else t.el=t.start:(t.start+n.ahead-1)
            tmp.el=matrix(marxss.model[["fixed"]][[el]][bad.names,1,t.el],length(bad.names),t.start+n.ahead-1)
            rownames(tmp.el)=bad.names
            newdata[[el]]=rbind(newdata[[el]],tmp.el)
          }
        }      
        #makesure the columns are in the same order as el.names
        newdata[[el]]=newdata[[el]][el.names,,drop=FALSE]      
        if(any(is.na(newdata[[el]]))){
          stop(paste("Stopped in predict_marxss(): There cannot be any NAs in the ",el," matrix.\n",sep=""),call.=FALSE)
        }
      }
    }
    if(("data" %in% names(newdata))){
      el="data"
      Y.names=attr(marxss.model,"Y.names")
      el.names=Y.names
      if(!is.matrix(newdata[[el]]))
        stop(c("Stopped in predict_marxss(): ", el," in newdata must be a matrix with time across the columns and n.ahead columns."),call.=FALSE)
      
      names.list=rownames(newdata[[el]])
      if(!all(el.names %in% names.list)){
        bad.names = el.names[!(el.names %in% names.list)]
        stop(c("Stopped in predict_marxss(): Some of the ",el," rows are missing: ",paste(bad.names,collapse=", ")),call.=FALSE)
      }
      #find any row names in the matrix that match the rownames in the el
      el.match=match(el.names,names.list)
      #subset those columns in the dataframe that match the el names
      el.in.matrix=el.names[!is.na(el.match)]
      newdata[[el]]=newdata[[el]][el.in.matrix,]
      #makesure the columns are in the same order as el.names
      newdata[[el]]=newdata[[el]][el.names,] 
    }
  }
  
  
  #now do some error checking
  if(!all(names(newdata) %in% allow.in)) stop(paste("Stopped in predict_marxss(): only allowed inputs for form ",form," are ",allow.in,".",sep=""),call.=FALSE)
  #Check that the dims of inputs are the same (or = 1) and fix c or d that have length 1 to have same dim 2 as other inputs
  thedims=unlist(lapply(newdata,function(x){dim(x)[2]}))
  #only consider those with length!=1
  thedims=thedims[thedims!=1]
  #n.ahead will have some value, since n.ahead is 1 by default in predict_marssMLE which called this function.
  #make sure n.ahead and any dim 2 of data, c or d !=1 are the same
  thedims=c(n.ahead,thedims)
  if(!all( abs(thedims - mean(thedims)) < .Machine$double.eps )){
    stop(paste("Stopped in predict_marxss(): in newdata, the 2nd dimension of all inputs (and n.ahead if passed in) must be equal (if not =1).",sep=""))
  }
  
  #set data to be all missing if it wasn't passed in
  if(is.null(newdata[["data"]])){
    newdata[["data"]] = matrix(as.numeric(NA), marxss.dims[["y"]][1], n.ahead)
    rownames(newdata[["data"]]) = attr(marxss.model,"Y.names")
  }else{ #it was passed in
    if(dim(newdata[["data"]])[2]!=n.ahead)
      stop(c("Stopped in predict_marxss(): data in newdata must be a matrix with time across the columns and n.ahead columns."),call.=FALSE)
  }
  
  #Now newdata ready as a list with elements c and d that correspond to model c and d as passed into MARSS()
  
  #### make a list of time-varying parameters
  param.names = attr(marxss.model,"par.names")
  time.varying.fixed = c()
  time.varying.free = c()
  for( elem in param.names ) {
    if( dim(marxss.model[["free"]][[elem]])[3] == 1 ){
      time.varying.free[elem] = FALSE  #not time-varying
    }else{ time.varying.free[elem] = TRUE }
    if( dim(marxss.model[["fixed"]][[elem]])[3] == 1 ){
      time.varying.fixed[elem] = FALSE  #not time-varying
    }else{ time.varying.fixed[elem] = TRUE }
  }
  time.varying = time.varying.fixed | time.varying.free
  
  if(any(time.varying)){
    if((t.start+n.ahead-1)>TT)
      stop(paste("Stopped in predict_marxss(): ",paste(param.names[time.varying], collapse=", ")," are time-varying.\nIn this case, you cannot forecast past the end of the time series\n(t.start+n.ahead must be < length of original data).\n",sep=""),call.=FALSE)
    param.t=t.start:(t.start+n.ahead-1)
  }
  #Now we need to construct a marssMODEL object (form=marxss) for predicting
  pred.marxss=marxss.model #copy the elements and attributes that won't change
  pred.model.dims=marxss.dims
  pred.model.dims[["data"]][2]=n.ahead
  pred.model.dims[["x"]][2]=n.ahead
  pred.model.dims[["y"]][2]=n.ahead
  for( el in param.names[!(param.names %in% c("c","d"))] ){
    if(!time.varying.fixed[el]){
      pred.marxss[["fixed"]][[el]]=marxss.model[["fixed"]][[el]]
    }else{ #is time.varying
      pred.marxss[["fixed"]][[el]]=marxss.model[["fixed"]][[el]][,,param.t,drop=FALSE]
      pred.model.dims[[el]][3]=n.ahead
    }
    if(!time.varying.free[el]){
      pred.marxss[["free"]][[el]]=marxss.model[["free"]][[el]]
    }else{ #is time.varying
      pred.marxss[["free"]][[el]]=marxss.model[["free"]][[el]][,,param.t,drop=FALSE]
      pred.model.dims[[el]][3]=n.ahead
    }  }
  
  for(el in c("c","d")){
    #this forces 2nd dim of matrix to be n.ahead
    pred.model.dims[[el]][3]=n.ahead
    #now replace any c or d row passed in with the new values
    pred.marxss[["fixed"]][[el]] = array(newdata[[el]],dim=c(dim(newdata[[el]])[1],1,n.ahead))
    rownames(pred.marxss[["fixed"]][[el]]) = rownames(newdata[[el]])
    #by definition free c and d is not estimated
    pred.marxss[["free"]][[el]] = marxss.model[["free"]][[el]]
  }
  
  #Set the initial conditons
  #Set these to the estimated distribution of x at t.start-1 conditioned on all the data (original)
  pred.marxss[["tinitx"]]=0
  kf=MARSSkf(x) #x is the original marssMLE obj passed in
  if(t.start==1){
    pred.marxss$fixed$x0=array(kf$x0T, dim=marxss.dims$x0)
    pred.marxss$fixed$V0=array(vec(kf$V0T), dim=c(marxss.dims$V0[1]*marxss.dims$V0[2],1,1))
  }else{
    pred.marxss$fixed$x0=array(kf$xtT[,t.start-1], dim=marxss.dims$x0)
    pred.marxss$fixed$V0=array(vec(kf$VtT[,,t=(t.start-1)]), dim=c(marxss.dims$V0[1]*marxss.dims$V0[2],1,1))
  }
  pred.marxss$free$x0=array(0,dim=c(marxss.dims$x0[1],0,1))
  pred.marxss$free$V0=array(0,dim=c(marxss.dims$V0[1]*marxss.dims$V0[2],0,1))
  
  pred.marxss[["data"]]=newdata[["data"]]
  attr(pred.marxss, "model.dims")=pred.model.dims #the marxss model dims
  
  ##Next step is to create the marssMLE object ready for predict.  This is done by changing x

  #These are not estimated when x_(tstart-1) and V_(t.start-1) are fixed at expected values
  for(val in c("par","start")){
    for(el in c("x0","V0")){
      if(!is.null(x[[val]][[el]])){
        x[[val]][[el]]=matrix(0,0,1)
      }
    }
  }
  #These don't have meaning when x_(tstart-1) and V_(t.start-1) are fixed at expected values
  for(val in c("par.se","par.bias","par.upCI","par.lowCI")){
    for(el in c("x0","V0")){
      if(!is.null(x[[val]][[el]])){
        x[[val]][[el]]=NULL
      }
    }
  }
  
  
  #now pred.marxss is a marssMODEL object (form=marxss) to be used for prediction
  x[["model"]]=pred.marxss
  #convert pred.marxss to a marss object and put in $marss
  x[["marss"]]=marxss_to_marss(pred.marxss)
  #only pass in marssMODEL since x par has been corrected above
  
  #x is now a marssMLE object with its marss and model changed to reflect newdata
  # with y, c and d potentially changed and t.start and t.end changed and TT changed
  # par$x0 and $V0 changed to reflect t.start
 
  #return this object to predict.marssMLE; only MLEobj is needed
  #but newdata list returned for debugging
  return(list(MLEobj=x, newdata=newdata))
}

describe_marxss = function(MODELobj){ describe_marss(MODELobj) }
#describe_marss works generally with marxss form models (of which marss is one type)
#describe_marss is in the file describe_marssMODEL.R


########################################################################
# is.marssMODEL_marxss function
# Check that the marxss object has all the parts it needs
# fixed, free, and par.names
# and that these have the proper size and form
# m is pulled from fixed$x0
########################################################################
is.marssMODEL_marxss <- function(MODELobj, method="kem"){
  
  msg=NULL

  if( !("marxss" %in% attr(MODELobj, "form") ) ){ 
    msg = c(msg, "Form attribute of the model object does not include marxss.\n")
  }
  ## Set up par.names that should be marxss model
  en = c("Z","A","R","B","U","Q","x0","V0","D","C","d","c","G","H","L")
  
  #Check that par.names has these and only these names
  par.names=attr(MODELobj, "par.names")
  if( !all(en %in% par.names ) ){ 
    msg = c(msg, "Element ", en[!(en %in% par.names)], " is missing from the par.names attribute of the model object.\n")
  }
  if( !all( par.names %in% en ) ) { 
    msg = c(msg, "Only ", en, "should be in the par.names attribute of the model object.\n")
  }
  if(!is.null(msg)){  #rest of the tests won't work so stop now
    return(msg)
  }
  
  ###########################
  # Check model.dims attribute is correct
  ###########################
  n = dim(MODELobj$data)[1]
  TT = dim(MODELobj$data)[2]
  m = dim(MODELobj$fixed$x0)[1]
  c1=dim(MODELobj$fixed$c)[1]
  d1=dim(MODELobj$fixed$d)[1]
  g1=dim(MODELobj$fixed$G)[1]/m
  h1=dim(MODELobj$fixed$H)[1]/n
  l1=dim(MODELobj$fixed$L)[1]/m
  
  en = c("Z","A","R","B","U","Q","x0","V0","D","C","d","c", "G", "H", "L", "data", "x", "y", "w", "v")
  correct.dim1 = c(Z=n,A=n, R=h1, B=m, U=m, Q=g1, x0=m, V0=l1, D=n, C=m, c=c1, d=d1, G=m, H=n, L=m, data=n, x=m, y=n, w=m, v=n)
  correct.dim2 = c(Z=m,A=1, R=h1, B=m, U=1, Q=g1, x0=1, V0=l1, D=d1, C=c1, c=TT, d=TT, G=g1, H=h1, L=l1, data=TT, x=TT, y=TT, w=TT, v=TT)
  model.dims=attr(MODELobj, "model.dims")
  for (elem in en) {
    ## Check for problems in the fixed/free pairs. Problems show up as TRUE 
    dim.flag1 = dim.flag2 = FALSE
    
    # check dim
    dim.flag1 = c(dim.flag1,!( model.dims[[elem]][1] == correct.dim1[[elem]] )) 
    dim.flag2 = c(dim.flag2,!( model.dims[[elem]][2] == correct.dim2[[elem]] ))
  }  
  if (any(c(dim.flag1, dim.flag2))) {  #There's a problem
    if(any(dim.flag1)) {
      msg = c(msg, paste("Dim 1 of ", en[dim.flag1], "is incorrect. Dims should be ", correct.dim1[dim.flag1], ", for a marss model.\n"))
    }
    if(any(dim.flag2)) {
      msg = c(msg, paste("Dim 2 of ", en[dim.flag2], "is incorrect. Dims should be ", correct.dim2[dim.flag2], ", for a marss model.\n"))
    }
    msg=c("\nErrors were caught in is.marssMODEL_marxss()\n", msg)
    return(msg)
  }  
  
  fixed=MODELobj$fixed; free=MODELobj$free
  
  ###########################
  # Check that x0, V0 and L are not time-varying
  ###########################
  en = c("x0", "V0", "L")
  time.var = NULL
  for (elem in en) {
    time.var.flag = FALSE
    time.var.flag = dim(fixed[[elem]])[3]!=1 | dim(free[[elem]])[3]!=1
    time.var <- c(time.var, time.var.flag)
  }
  if(any(time.var)) {  #There's a problem
    msg = c(msg, paste(en[time.var], "cannot be time-varying.  3rd dim of fixed and free must equal 1.\n"))
    msg=c("\nErrors were caught in is.marssMODEL_marxss()\n", msg)
    return(msg)
  }
  
  ###########################
  # Check that none of the var-cov matrices have negative values on the diagonal
  # and that there are no f+Dq elements only f+0q or 0+Dq
  # and D must be a design matrix, so no beta_1*q1 + beta_2*q2 elements
  ###########################
  en = c("R", "Q", "V0")
  neg = bad.var = not.design = NULL
  for (elem in en) {
    neg.flag = bad.var.flag = not.design.flag = FALSE
    for(i in 1:max(dim(free[[elem]])[3],dim(fixed[[elem]])[3])){
      if(dim(fixed[[elem]])[3]==1){i1=1}else{i1=i}
      if(dim(free[[elem]])[3]==1){i2=1}else{i2=i}
      if(is.fixed(free[[elem]][,,min(i,dim(free[[elem]])[3]),drop=FALSE])){ #this works on 3D mats
        zero.free.rows = matrix(TRUE,correct.dim1[[elem]]*correct.dim2[[elem]],1)
      }else{
        zero.free.rows=apply(free[[elem]][,,i2,drop=FALSE]==0,1,all) #works on 3D mat
        #the requirement is that each estimated element (in p) appears only in one place in the varcov mat, but fixed rows (0 rows) are ok
        not.design.flag = !is.design(free[[elem]][,,i2,drop=FALSE],strict=FALSE,zero.rows.ok=TRUE,zero.cols.ok=TRUE) #works on 3D if dim3=1
      }
      zero.fixed.rows=apply(fixed[[elem]][,,i1,drop=FALSE]==0,1,all) #works on 3D
      fixed.mat = unvec(fixed[[elem]][,,i1],dim=c(correct.dim1[[elem]],correct.dim2[[elem]]))
      if( any(!zero.fixed.rows & !zero.free.rows) ){
        bad.var.flag = TRUE   #no f+Dq rows
      }
      if(any(takediag(fixed.mat)<0,na.rm=TRUE)) neg.flag=TRUE      #no negative diagonals
    } #end the for loop over time
    not.design = c(not.design, not.design.flag)
    neg = c(neg, neg.flag)
    bad.var = c(bad.var, bad.var.flag)
  } #enf the for loop over elem
  if(any(neg)) {
    msg = c(msg, paste("Negative values are on the diagonal of ", en[neg], ". Neg values are illegal on the diag of a var-cov matrix.\n", sep=""))
  }
  if(any(bad.var)) {
    msg = c(msg, paste("Fixed and estimated values are combined in some elements of ", en[bad.var], ". This is not allowed.\n", sep=""))
  }
  if(any(not.design)) {
    msg = c(msg, paste("The D matrices of ", en[not.design], " must be design matrices.\n", sep=""))
  }
  
  ###########################
  # Check that V0, Q and R matrices are symmetric and positive-definite
  ###########################
  en = c("R", "Q", "V0")
  pos = symm = NULL
  for (elem in en) {
    varcov.flag = TRUE; varcov.msg=""
    var.dim = c(correct.dim1[[elem]],correct.dim2[[elem]])
    for(i in 1:model.dims[[elem]][3]){
      if(dim(fixed[[elem]])[3]==1){i1=1}else{i1=i}
      if(dim(free[[elem]])[3]==1){i2=1}else{i2=i}
      #works on 3D if dim3=1
      par.as.list = fixed.free.to.formula(fixed[[elem]][,,i1,drop=FALSE],free[[elem]][,,i2,drop=FALSE],var.dim) #coverts the fixed,free pair to a list matrix
      tmp=is.validvarcov(par.as.list, method=method)
      varcov.flag=varcov.flag & tmp$ok
      if(!tmp$ok) varcov.msg = c(varcov.msg, paste(" ", tmp$error, "at t=", i, "\n",sep=""))
      
      if(!varcov.flag) msg = c(msg, paste("The variance-covariance matrix ", elem, " is not properly constrained.\n", sep=""), varcov.msg)
    } #end for loop over time

  } #end for loop over elements

  ###########################
  # Check that crossprod(G), crossprod(H), crossprod(L) are invertible
  ###########################
  en = c("G", "H", "L")
  pos = symm = NULL
  for (elem in en) {
    varcov.flag = TRUE; varcov.msg=""
    var.dim = c(correct.dim1[[elem]],correct.dim2[[elem]])
    for(i in 1:model.dims[[elem]][3]){
      if(dim(fixed[[elem]])[3]==1){i1=1}else{i1=i}
      if(dim(free[[elem]])[3]==1){i2=1}else{i2=i}
      #works on 3D if dim3=1
      #since G, H, and L are numeric, par.as.list will be a numeric matrix not list
      par.as.list = fixed.free.to.formula(fixed[[elem]][,,i1,drop=FALSE],free[[elem]][,,i2,drop=FALSE],var.dim) #coverts the fixed,free pair to a list matrix
      #this requirement is mention in 4.4 in EM Derivation
      #simple test for invertibility via condition number
      condition.limit=1E10
      tmp=kappa(crossprod(as.numeric(par.as.list))) < condition.limit #TRUE is good
      varcov.flag=varcov.flag & tmp
      if(!tmp) varcov.msg = c(varcov.msg, paste(" ", tmp$error, "at t=", i, "\n",sep=""))
      
      if(!varcov.flag) msg = c(msg, paste("The matrix t(", elem, ")%*%", elem," must be invertible.\n", sep=""), varcov.msg)
    } #end for loop over time
    
  } #end for loop over elements
  
if(length(msg) == 0){ return(NULL)
}else {
  msg=c("\nErrors were caught in is.marssMODEL_marxss()\n", msg)
  return(msg)
}
}

