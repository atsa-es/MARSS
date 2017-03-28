MARSSapplynames=function(MLEobj){
## Helper function to put names on the elements in a marssMLE object
  if(!(class(MLEobj) %in% c("marssMLE")))
     stop("Stopped in MARSSapplynames() because this function is for marssMLE objects only.\n", call.=FALSE)

  modelObj=MLEobj[["marss"]]
  par.names=attr(modelObj,"par.names")
  X.names=attr(modelObj,"X.names")
  Y.names=attr(modelObj,"Y.names")
  
  #The par element has rownames that come from the column names of the free matrix   
  for(elem in par.names){
    if(!is.null(MLEobj[["par"]][[elem]]) & is.null(rownames(MLEobj[["par"]][[elem]]))) rownames(MLEobj[["par"]][[elem]]) = colnames(modelObj$free[[elem]])
  }
  
  rownames(MLEobj[["model"]][["data"]])=Y.names
  rownames(MLEobj[["marss"]][["data"]])=Y.names
  if(!is.null(MLEobj[["Ey"]][["ytT"]])) rownames(MLEobj[["Ey"]][["ytT"]]) =  Y.names
  if(!is.null(MLEobj[["kf"]][["xtT"]])) rownames(MLEobj[["kf"]][["xtT"]]) =  X.names
  if(!is.null(MLEobj[["kf"]][["xtt1"]])) rownames(MLEobj[["kf"]][["xtt1"]]) =  X.names
  if(!is.null(MLEobj[["states.se"]])) rownames(MLEobj[["states.se"]]) =  X.names
  if(!is.null(MLEobj[["states"]])) rownames(MLEobj[["states"]]) =  X.names
  
  return(MLEobj)
  
}
