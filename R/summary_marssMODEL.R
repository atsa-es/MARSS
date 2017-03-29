###############################################################################################################################################
#  Summary method for class marssMODEL
###############################################################################################################################################

summary.marssMODEL <- function (object, ..., silent = FALSE) 
{

   n = dim(object$data)[1]; m = dim(object$fixed$x0)[1]
   cat(paste("Model Structure is\n","m: ",m," state process(es)\n","n: ",n," observation time series\n",sep=""))

   rpt.list = list()
   en = attr(object,"par.names") 
   dim.tmp = attr(object, "model.dims") 

   xnames = attr(object,"X.names")
	 if(is.null(xnames)) xnames = paste("X",1:m,sep="")
   ynames = attr(object,"Y.names")
   if(is.null(ynames)) ynames = paste("Y",1:m,sep="")
   
  for (elem in en) {
    #list matrix version of the model    
    Tmax=dim.tmp[[elem]][3]
    for(t in 1:Tmax){
      tmp = fixed.free.to.formula( sub3D(object$fixed[[elem]],t=t),sub3D(object$free[[elem]],t=t),dim.tmp[[elem]][1:2] )
      if(Tmax==1) rpt.list[[elem]] = tmp
      if(t==1 & Tmax>1) rpt.list[[elem]]=array(list(),dim=c(dim(tmp),Tmax))
      if(Tmax>1) rpt.list[[elem]][,,t] = tmp
    } # for t in Tmax
	  if(elem %in% c("B","Q","Z","V0")) colnames(rpt.list[[elem]])=xnames
	  if(elem %in% c("x0","U")) colnames(rpt.list[[elem]])="X"
	  if(elem %in% c("R")) colnames(rpt.list[[elem]])=ynames
	  if(elem %in% c("A")) colnames(rpt.list[[elem]])="Y"
	  if(elem %in% c("B","U","Q","x0","V0")) rownames(rpt.list[[elem]])=xnames
	  if(elem %in% c("Z","A","R")) rownames(rpt.list[[elem]])=ynames
  } #for elem
  rpt.list$tinitx = object$tinitx
   
  if(!silent) print(rpt.list, quote=FALSE)
  invisible(rpt.list)
          
  }