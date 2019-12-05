############################################################################################################################
#   checkModelList()
#   This checks user model list passed to MARSS().
#   The main purpose is to make sure that MARSS.form functions will work not to make sure model is valid
#   Dim checks on matrices are carried out in the MARSS.form functions that translate a model list to marssMODEL object
#   No error checking is done on controls and inits besides checking that it is present (NULL is ok); 
#   is.marssMLE() will error-check controls and inits
##########################################################################################################################
checkModelList = function( model, defaults, this.form.allows)
{
## First deal with case where model is not passed in all
if(is.null(model)) model=defaults

if (!is.list(model)){
  msg=" model must be passed in as a list.\n"
      cat("\n","Errors were caught in checkModelList \n", msg, sep="")
      stop("Stopped in checkModelList() due to specification problem(s).\n", call.=FALSE)
}
model.elem=names(defaults) 
### If some elements are missing from the model list use the defaults    
passed.in = model.elem %in% names(model)
for(el in model.elem[!passed.in] ){
    model[[el]] = defaults[[el]]
}
for(el in model.elem[passed.in] ){
    if(is.null(model[[el]])) model[[el]] = defaults[[el]]
}
#Check model structures (b497)
if( !all(names(model) %in% model.elem ) ){ 
      bad.name = names(model)[!(names(model) %in% model.elem )]
      msg=paste(" Elements ", bad.name," not allowed in model list for this form.\n",sep="")
      cat("\n","Errors were caught in checkModelList \n", msg, sep="")
      stop("Stopped in checkModelList() due to specification problem(s).\n", call.=FALSE)
      }
if( !all(model.elem %in% names(model)) ){ 
  bad.name = model.elem[!(model.elem %in% names(model))]
  msg=paste(" Element ", bad.name, " is missing in the model list passed into MARSS().\n", sep="")
      cat("\n","Errors were caught in checkModelList \n", msg, sep="")
      stop("Stopped in checkModelList() due to specification problem(s).\n", call.=FALSE)
      }

#check that model list doesn't have any duplicate names
  if(any(duplicated(names(model)))){
    bad.name = names(model)[duplicated(names(model))]
      msg=paste(" The elements ", bad.name, " are duplicated in the model list passed into MARSS().\n", sep="")
      cat("\n","Errors were caught in checkModelList \n", msg, sep="")
      stop("Stopped in checkModelList() due to specification problem(s).\n", call.=FALSE)
   }

#Series of checks on the model specification
problem = FALSE
msg=NULL
#check model structures only have allowed cases
  for (el in model.elem) {
    bad.str = FALSE
    #if length=1, then it must be a character or numeric string and that string must be in allowed. vectors length>1 are not allowed

    if(!is.factor(model[[el]]) && !is.array(model[[el]])) {
        if(length(model[[el]])!=1) bad.str=TRUE
        if(!bad.str){
          testit = try(model[[el]] %in% this.form.allows[[el]])
          if(inherits(testit,"try-error") ){ bad.str=TRUE
          }else{ if(!testit ) bad.str=TRUE }
        }
    }
    if(bad.str) {
      problem=TRUE
      msg = c(msg, paste(" The model value for ", el, " is not allowed. Check ?MARSS.form \n", sep=""))
      }
    #if factor, must be allowed to be factor
    if(is.factor(model[[el]]) && !(el %in% this.form.allows$factors)) {
      problem=TRUE
      msg = c(msg, paste(" model$",el," is not allowed to be a factor.\n", sep=""))
      }
    #if matrix must be allowed to be matrix
    if(is.array(model[[el]]) && !(el %in% this.form.allows$matrices)){
      problem=TRUE
      msg = c(msg, paste(" model$",el," is not allowed to be a matrix.\n", sep=""))
      }     
    #if matrix then no NAs if character or list; this would be caught in is.marssMODEL but would be hard for user to understand problem
    if(is.array(model[[el]]) && (el %in% this.form.allows$matrices)){
      if( any(is.na(model[[el]])) || any(unlist(lapply(model[[el]],is.infinite))) ){
        problem=TRUE
        msg = c(msg, paste(" model$",el," is a matrix. No NAs or Infs allowed in this case.\n", sep=""))
        }
      if( is.list(model[[el]]) && any(sapply(model[[el]],length)!=1) ){
        problem=TRUE
        msg = c(msg, paste(" model$",el," is a list matrix. Each element must be length 1.\n", sep=""))
        }
      } # is matrix  
    } # end for (el in model.elem)

#If there are errors, then don't proceed with the rest of the checks
  if(problem)  {
          cat("\n","Errors were caught in checkModelList \n", msg, sep="")
          stop("Stopped in checkModelList() due to specification problem(s).\n", call.=FALSE)
        }
        
  #return the updated model list with missing elements filled in with defaults
  model
}
