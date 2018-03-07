###############################################################################################################################################
#  model.frame method for class marssMLE and class marssMODEL 
#  returns a data.frame if model is in certain forms
##############################################################################################################################################
model.frame.marssMODEL <- function (formula, ...) {
  model=formula
  allowed.forms = c("dlm")
  model.form = attr(model,"form")[1]
  if(!(model.form %in% allowed.forms)) stop(paste("model.frame.marssMODEL doesn't know how to deal with marssMODEL objects of form ",model.form,".",sep=""))
  ret = data.frame(t(model$data))
  
  if(model.form=="dfa"){
    if(is.null(MLEobj$call)) stop("model.frame.marssMLE needs to have the call element in the marssMLE object.  This is used for the covariates.")
    covariates = MLEobj$call$covariates
    if(!is.null(covariates)){
      if(is.null(rownames(covariates))) rownames(covariates)=rownames(model$fixed$d)
      ret = cbind(ret, t(covariates))
    }
  }
  ret
}  #end of model.frame.marssMODEL

model.frame.marssMLE <- function (formula, ...) {
  model.frame.marssMODEL(formula$model)
}