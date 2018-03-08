###################################################################################
# Helper functions to create and work with a DLM model sensu Zuur
# x(t)=x(t-1) + w(t), W~MVN(0,1)
# y(t)=Z x(t) + A(t) + D(t) d(t) + v(t), V~MVN(0,R)
# x(t0) = x0 + l, L ~ MVN(0,5)

# MARSS.dfa: The conversion functions have 2 parts
# Part 1 Set up the DLM model in MARSS.marxss form
# Part 2 Call MARSS.marxss to finish the set-up and checking

# Functions that are called by the generic functions print.marssMLE, residuals.marssMLE,
# predict.marssMLE, coef.marssMLE, MARSSinits.marssMLE
# print_dfa, residuals_dfa, predict_dfa, coef_dfa, MARSSinits_dfa
###################################################################################
MARSS.dfa=function(MARSS.call){
  # MARSS(data, model=list(), covariates=NULL, z.score=TRUE, demean=TRUE, control=list())
  # model.defaults =list(A="zero", R="diagonal and equal", D="zero", x0="zero", V0=diag(5,1), tinitx=0, diffuse=FALSE, m=1)
  
  #load needed package globals
  common.allowed.in.MARSS.call=get("common.allowed.in.MARSS.call", envir=pkg_globals)
  
  
  #Part 1 Set up defaults and check that what the user passed in is allowed
  # 1 Check for form dependent user inputs for method and reset defaults for inits, and control if desired
  # 2 Specify the text shortcuts and whether factors or matrices can be passed in
  #   The names in the allowed list do not need to be A, B, Q .... as used in the marss form
  #   Other names can be used if you want the user to use those names; then in the MARSS.form function,
  #   you convert the user passed in names into the marss form with the A, B, Q, R, ... names
  #   checkModelList() will check what the user passes in against these allowed values, so
  #   so you need to make sure each name in model.defaults has a model.allowed value here 
  
  #Check that no args were passed into MARSS that are not allowed
  dfa.allowed.in.MARSS.call=c("model","z.score","demean","covariates")
  allowed.in.call=c(dfa.allowed.in.MARSS.call,common.allowed.in.MARSS.call)
  if(any(!(names(MARSS.call) %in% allowed.in.call))){
    bad.names=names(MARSS.call)[!(names(MARSS.call) %in% allowed.in.call)]
    msg=paste("Argument ", paste(bad.names, collapse=", "),"  not allowed MARSS call for form ", MARSS.call$form, ". See ?MARSS.dfa\n", sep="")
    cat("\n","Errors were caught in MARSS.dfa \n", msg, sep="")
    stop("Stopped in MARSS.dfa() due to problem(s) with model specification.\n", call.=FALSE)
  }
  
  #Set up some defaults
  if(is.null(MARSS.call[["z.score"]])) MARSS.call$z.score=TRUE
  if(is.null(MARSS.call[["demean"]])) MARSS.call$demean=TRUE
  
  #Start error checking
  problem=FALSE
  msg=c()
  #check that data and covariates elements are matrix or vector, no dataframes, and is numeric
  for(el in c("data","covariates"[!is.null(MARSS.call[["covariates"]])])){
    if( !(is.matrix(MARSS.call[[el]]) || is.vector(MARSS.call[[el]])) ) {
      problem=TRUE
      msg = c(msg, paste(el," must be a matrix or vector (not a data frame or 3D array).\n"))
    }else{ 
      if(is.vector(MARSS.call[[el]])) MARSS.call[[el]]=matrix(MARSS.call[[el]],1)
    }
  }
  if(!is.null(MARSS.call[["covariates"]])){
    if(dim(MARSS.call[["data"]])[2]!= dim(MARSS.call[["covariates"]])[2]){
      problem=TRUE
      msg = c(msg, "data and covariates must have the same number of time steps.\n")
    }
  }
  if(!is.null(MARSS.call[["covariates"]])){
    if(any(is.na(MARSS.call[["covariates"]]))){
      problem=TRUE
      msg = c(msg, "covariates cannot have any missing values in a standard DFA.\n See User Guide section on DFA for alternate approaches when covariates have missing values.\n")
    }
  }
  if(!is.null(MARSS.call[["model"]][["m"]])){
    if(length(MARSS.call[["model"]][["m"]])!=1){
      problem=TRUE
      msg = c(msg, "model$m must be an integer between 1 and n.\n")
    }else{
      if( !is.numeric(MARSS.call[["model"]][["m"]]) | !is.wholenumber(MARSS.call$model$m) ){
        problem=TRUE
        msg = c(msg, "model$m must be an integer between 1 and n.\n")
      }else{
        if(MARSS.call[["model"]][["m"]] > dim(MARSS.call[["data"]])[1]){
          problem=TRUE
          msg = c(msg, "model$m must be an integer between 1 and n.\n")
        }
      }
    } }
  
  if(problem){
    cat("\n","Errors were caught in MARSS.dfa \n", msg, sep="")
    stop("Stopped in MARSS.dfa() due to specification problem(s).\n", call.=FALSE)
  }
  
  n=dim(MARSS.call[["data"]])[1]
  model.allowed = list(
    #if it is a length 1 vector then the value must be one of these.  All elements in your model list must be here
    A=c("unequal","zero"),
    R=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    D=c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    x0=c("unconstrained", "unequal", "zero"),
    B=c("identity","diagonal and equal","diagonal and unequal"),
    Q=c("identity","diagonal and equal","diagonal and unequal"),
    V0=c("identity", "zero"),
    tinitx=c(0,1),
    diffuse=c(TRUE,FALSE),
    m=1:n,
    #This line says what is allowed to be a matrix
    matrices = c("Z","A","R","D","x0","V0")
  )
  
  if(is.null(MARSS.call[["model"]][["m"]])) m=1 else m=MARSS.call[["model"]][["m"]]
  if(!is.null(MARSS.call[["covariates"]])) D="unconstrained" else D="zero"
  if(is.null(MARSS.call[["model"]][["Z"]])){ #Set up default Z
    Z = matrix(list(), nrow=n, ncol=m)
    # insert row (i) & col (j) indices
    for(i in seq(n)) {Z[i,] = paste(i, seq(m), sep="")}
    # set correct i,j values in Z to numeric 0
    if(m>1) for(i in 1:(m-1)){Z[i,(i+1):m] = 0}
  }
  #defaults for any missing model list elements
  model.defaults =list(
    A="zero", 
    R="diagonal and equal", 
    D=D, 
    B="identity",
    Q="identity",
    x0="zero", 
    V0=diag(5,m), 
    tinitx=0, 
    diffuse=FALSE, 
    m=1,
    Z=Z)
  
  #This checks that what user passed in model list can be interpreted and converted to form marss
  #if no errors, it updates the model list by filling in missing elements with the defaults
  MARSS.call$model=checkModelList( MARSS.call[["model"]], model.defaults, model.allowed )
  model=MARSS.call[["model"]]
  
  if(!(MARSS.call$z.score %in% c(TRUE,FALSE)))
    stop("Stopped in MARSS.dlm(): z.score must be TRUE/FALSE.\n", call.=FALSE)
  if(!(MARSS.call$demean %in% c(TRUE,FALSE)))
    stop("Stopped in MARSS.dlm(): demean must be TRUE/FALSE.\n", call.=FALSE)
  
  #Set up U; always 0 for dfa
  U=matrix(0,m,1)
  
  #Set up D and d
  if(is.null(MARSS.call[["covariates"]])) d=matrix(0,1,1) else d=MARSS.call$covariates
  
  #Set up Q  & B; always fixed for dfa
  Q=diag(1,m); B=diag(1,m)
  
  # set up list of model components for a marxss model
  dfa.model = list(Z=Z, A=model$A, D=model$D, d=d, R=model$R, B=B, U=U, Q=Q, x0=model$x0, V0=model$V0, tinitx=model$tinitx)
  
  dat=MARSS.call[["data"]]
  if(MARSS.call$demean){
    y.bar = apply(dat, 1, mean, na.rm=TRUE)
    dat = (dat - y.bar)
  }
  if(MARSS.call[["z.score"]]){
    Sigma = sqrt(apply(dat, 1, var, na.rm=TRUE))
    dat = dat * (1/Sigma)
  }
  
  MARSS.call = list(data=dat, inits=MARSS.call$inits, control=MARSS.call$control, method=MARSS.call$method, form="dlm", silent=MARSS.call$silent, fit=MARSS.call$fit, fun.kf=MARSS.call$fun.kf)
  
  #dfa is a type of marxss model, so use MARSS.marxss to test it and set up the marss object
  tmp = MARSS.marxss(list(data=dat,model=dfa.model,method=MARSS.call$method,silent=MARSS.call$silent))
  #marss is the name for the form=marss model object that MARSS.form functions return
  #need to add "dfa" to attribute form
  marxss_object=tmp$model
  attr(marxss_object, "form")=c("dfa", "marxss")
  MARSS.call$model = marxss_object
  MARSS.call$marss = tmp$marss
  
  ## Return MARSS inputs as list
  MARSS.call
}
#This works since dfa just creates a marxss object so x$model is marssMODEL form=c("dfa", "marxss")
#probably want to customize for dfa later
print_dfa = function(x){ return(print_marxss(x)) }
coef_dfa = function(x){ return(coef_marxss(x)) }
MARSSinits_dfa = function(MLEobj, inits){ return(MARSSinits_marxss(MLEobj, inits)) }
predict_dfa = function(x, newdata, n.ahead, t.start){ predict_marxss(x, newdata, n.ahead, t.start) }
describe_dfa = function(MODELobj){ describe_marss(MODELobj) }
is.marssMODEL_dfa = function(MODELobj, method="kem"){ is.marssMODEL_marxss(MODELobj, method=method) }
augment_dfa = function(x, type.predict, interval, conf.level, extra){
  if(is.null(extra[["rotate"]])){ rotate=FALSE }else{ rotate=extra[["rotate"]] }
  
  ret = augment_marxss(x, type.predict=type.predict, interval=interval, conf.level=conf.level)
  
  if(type.predict=="states")
    stop("Augment does not return the estimated states (loadings).  See ?augment.marssMLE . Use tidy().\n")
  
  if(rotate){
    coefs = coef(MLEobj, type="matrix")
    ZZ = coefs[["Z"]]
    DD = coefs[["DD"]]
    dd = coefs[["dd"]]
    AA = coefs[["A"]]
    H_inv <- varimax(ZZ)[["rotmat"]] # inv of the rotation matrix
    ## check for covars
    if(!is.null(covariates)) {
      DD = coef(MLEobj, type="matrix")$D
      cov_eff = DD %*% covariates
    } else {
      cov_eff = matrix(0, nn, TT)
    }
    ret$.fitted = ZZ %*% H_inv %*% MLEobj$states + DD %*% dd + AA
    if(interval=="confidence"){
      alpha = 1-conf.level
      ret$.conf.low = qnorm(alpha/2)*ret$.se.fit + ret$.fitted
      ret$.conf.up = qnorm(1-alpha/2)*ret$.se.fit + ret$.fitted
    }
  }
  ret
}
