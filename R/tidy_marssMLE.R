###############################################################################################################################################
#  tidy method for class marssMLE
##############################################################################################################################################
tidy.marssMLE = function (x,  type = c("parameters", "states"), 
                          conf.int = TRUE, conf.level = 0.95, ...) 
{ 
  type = match.arg(type)
  alpha = 1-conf.level
  
  if(type=="parameters"){
  ests = coef(x, type="vector")
  if(length(ests)==0) stop("tidy.marssMLE: No estimated parameters in your fitted model.\n")

  model.has.cis = all(c("par.se", "par.lowCI", "par.upCI")%in%names(x))
  if(conf.int) rerun.MARSSparamCIs = ifelse(model.has.cis, FALSE, TRUE)
  if(!missing(...)){
    extras=list(...)
    if(!all(names(extras)%in%c("method", "hessian.fun", "nboot"))) stop("Unknown extra argument. See ?tidy.marssMLE for allowed arguments.\n")
  }
  if(conf.int & (!missing(...) | !missing(conf.level))){
    if(model.has.cis) warning("tidy.marssMLE: Your marssMLE object has par.se and CIs, but you have passed in arguments for calculating CIs.  MARSSparamCIs() will be re-run with these values.\n")
    rerun.MARSSparamCIs = TRUE
  }
  ret = data.frame(
    term = names(ests),
    estimate = ests
  )
  if( conf.int ){
    if(rerun.MARSSparamCIs) x = MARSSparamCIs(x, ...)
    ret = cbind(ret, 
            std.error = coef(x, type="vector", what="par.se"),
            conf.low = coef(x, type="vector", what="par.lowCI"),
            conf.high = coef(x, type="vector", what="par.upCI")
    )
  }
  rownames(ret)=NULL
  }
  if(type=="states"){
    model=x[["model"]]
    state.names = attr(model, "X.names")
    state.dims = attr(model, "model.dims")[["x"]]
    mm = state.dims[1]
    TT = state.dims[2]
    ret = data.frame(
      term=rep(state.names,each=TT), 
      t=rep(1:TT,mm), 
      estimate = vec(t(x[["states"]])),
      std.error = vec(t(x[["states.se"]]))
    )
    if( conf.int ){
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
