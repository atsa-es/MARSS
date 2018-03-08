###############################################################################################################################################
#  augment method for class marssMLE
#  returns fitted values, residuals, std err of residuals and std residuals
#  the base function is augment_marxss and is in MARSS_marxss.R
#  augment_dfa uses the optional rotate argument
##############################################################################################################################################
augment.marssMLE <- function (x, type.predict = c("observations", "states"),
                              interval = c("none", "confidence"), conf.level = 0.95, ...) {
  args=list(...)
  if(length(args)!=0){
  if(!all(names(args)%in%c("rotate"))) stop("Unknown extra argument. See ?augment.marssMLE for allowed arguments.\n")
  }
  
  type.predict = match.arg(type.predict)
  interval = match.arg(interval)

  form=attr(x[["model"]], "form")
  augment.fun = paste("augment_", form[1], sep="")
  tmp=try(exists(augment.fun, mode="function"),silent=TRUE)
  if(isTRUE(tmp)){
    ret=eval(call(augment.fun, x, type.predict = type.predict, interval = interval, conf.level = conf.level, extra=args))
  }else{ 
    ret=c(msg, paste("No augment_", form[1], " is available.\n", sep=""))
  }
  ret
}
