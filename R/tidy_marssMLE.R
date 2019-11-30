###############################################################################################################################################
#  tidy method for class marssMLE
##############################################################################################################################################
tidy.marssMLE = function (x,  type = c("parameters", "states", "observations", "xtT", "xtt1", "xtt", "ytT"), 
                          interval = TRUE, alpha = 0.05,
                          form=attr(x[["model"]], "form")[1], ...)
{ 
  ## Argument checking
  type = match.arg(type)
  if(!is.numeric(alpha) | alpha > 1 | alpha < 0) stop("tidy.marssMLE: alpha must be between 0 and 1.", call. = FALSE)
  if(!(interval%in%c(TRUE,FALSE))) stop("tidy.marssMLE: interval must be TRUE/FALSE", call. = FALSE)
  if(type=="states") type="xtT"
  if(type=="observations") type="ytT"
  ## End Argument checking
  
  extras=list()
  
  rerun.MARSSparamCIs = FALSE
  model.has.cis = all(c("par.se", "par.lowCI", "par.upCI")%in%names(x))
  if(interval & type=="parameters") rerun.MARSSparamCIs = ifelse(model.has.cis, FALSE, TRUE)
  if(!missing(...)){
    extras=list(...)
    if(!all(names(extras)%in%c("rotate", "method", "hessian.fun", "nboot"))) stop("Unknown extra argument. See ?tidy.marssMLE for allowed arguments.\n")
  }
  
  if(type=="parameters" & interval & (!missing(...) | !missing(alpha))){
    if(model.has.cis) warning("tidy.marssMLE: Your marssMLE object has par.se and CIs, but you have passed in arguments for calculating CIs.  MARSSparamCIs() will be re-run with these values.\n")
    rerun.MARSSparamCIs = TRUE
  }
  
  #set rotate
  if(!(form%in%c("marss","marxss","dfa")))  stop("tidy.marssMLE: Allowed forms are marss, marxss, and dfa.\n", call.=FALSE)
  rotate = FALSE
  if(form=="dfa" & "rotate"%in%names(extras)){
    rotate=extras[["rotate"]]
    if(!(rotate %in% c(TRUE, FALSE))) stop("tidy.marssMLE: rotate must be TRUE/FALSE. \n", call.=FALSE)
  }
  if(form!="dfa" & "rotate"%in%names(extras)) 
    if(rotate) stop("tidy.marssMLE: rotate only makes sense if form='dfa'.\n  Pass in form='dfa' if your model is a DFA model, but the form\n attribute is not set (because you set up your DFA model manually). \n", call.=FALSE)
  
  if(type=="parameters"){
    ests = coef(x, type="vector")
    if(length(ests)==0) stop("tidy.marssMLE: No estimated parameters in your fitted model.\n", call.=FALSE)
    if(form=="dfa" & rotate & length(x[["par"]][["Z"]])!=0){
      stop("tidy.marssMLE: You are requesting the parameters for a DFA \n and requested that the Z matrix be rotated. You need to do the rotation yourself.  See ?tidy.marssMLE for the code.\n", call.=FALSE)
    }else{
      ret = data.frame(
        term = names(ests),
        estimate = ests
      )
      if( interval ){
        if(rerun.MARSSparamCIs) x = MARSSparamCIs(x, alpha=alpha, ...)
        ret = cbind(ret, 
                    std.error = coef(x, type="vector", what="par.se"),
                    conf.low = coef(x, type="vector", what="par.lowCI"),
                    conf.high = coef(x, type="vector", what="par.upCI")
        )
      }
    }
    rownames(ret)=NULL
  }
  if(type%in%c("xtt1","xtT","xtt")){
    model=x[["model"]]
    state.names = attr(model, "X.names")
    state.dims = attr(model, "model.dims")[["x"]]
    mm = state.dims[1]
    TT = state.dims[2]
    kfss=MARSSkfss(x)
    states = kfss[[type]]
    vtype = str_replace(type,"x","V")
    states.se = apply(kfss[[vtype]],3,function(x){takediag(x)})
    states.se[states.se<0]=NA
    states.se=sqrt(states.se)
    if(m==1) states.se=matrix(states.se,1,TT)
    rownames(states.se)=attr(MLEobj$marss, "X.names")
    #if user specified rotate
    if(form=="dfa" & rotate & length(x[["par"]][["Z"]])!=0){
      Z.est = coef(x, type="matrix")[["Z"]]
      H = 1
      if(ncol(Z.est)>1){
        H = solve(varimax(Z.est)[["rotmat"]])
        states = H %*% states #rotated states
        TT = attr(x$model, "model.dims")[["data"]][2]
        states.var=MARSSkf(x)[[vtype]]
        for(t in 1:TT){
          states.se[,t] = sqrt(takediag(H%*%states.var[,,t]%*%t(H)))
        }
      }
    }
    ret = data.frame(
      term=rep(state.names,each=TT), 
      t=rep(1:TT,mm), 
      estimate = vec(t(states)),
      std.error = vec(t(states.se))
    )
    if( interval ){
      conf.low = qnorm(alpha/2)*ret$std.error + ret$estimate
      conf.up = qnorm(1-alpha/2)*ret$std.error + ret$estimate
      ret = cbind(ret, 
                  conf.low = conf.low,
                  conf.high = conf.up
      )
    }
    rownames(ret)=NULL
  }
  if(type%in%c("ytT")){
    model=x[["model"]]
    Y.names = attr(model, "Y.names")
    Y.dims = attr(model, "model.dims")[["y"]]
    nn = Y.dims[1]
    TT = Y.dims[2]
    aug <- broom::augment(x, type.predict="observations")

    ret = data.frame(
      term=rep(Y.names,each=TT), 
      t=rep(1:TT,nn), 
      estimate = aug$.fitted,
      std.error = aug$.se.fit
    )
    if( interval ){
      conf.low = qnorm(alpha/2)*ret$std.error + ret$estimate
      conf.up = qnorm(1-alpha/2)*ret$std.error + ret$estimate
      ret = cbind(ret, 
                  conf.low = conf.low,
                  conf.high = conf.up
      )
    }
    rownames(ret)=NULL
  }
  ret
}
