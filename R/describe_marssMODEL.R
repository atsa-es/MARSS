#This defines the describe.marssMODEL which calls in turn describe_form functions
describe.marssMODEL = function(x){
  if( !("marssMODEL" %in% class(x)) ) stop("Stopped in describe.marssMODEL(): x must be a marssMODEL object.\n",call.=FALSE)
  
  form=attr(x,"form")
  #First make sure specified equation form has a corresponding function to do the conversion to marssMODEL (form=marss) object
  describe.fun = paste("describe_",form[1],sep="")
  #must return a list with the model described
  tmp=try(exists(describe.fun,mode="function"),silent=TRUE)
  if(!isTRUE(tmp)){
    msg=paste("describe.marssMODEL: describe_", form[1], "() function to describe the marssMODEL (form=",form[1],") does not exist.\n",sep="")
    stop(msg, call.=FALSE)
  }
  constr.type = eval(call(describe.fun, x))
  return(constr.type)
}
