###############################################################################################################################################
#  Predict method for class marssMLE. 
##############################################################################################################################################
predict.marssMLE <- function (object, ..., n.ahead=1, t.start=NULL, newdata=list(), se.fit=TRUE, nboot=1000, param.gen="hessian", verbose=FALSE, prediction.intervals=TRUE) {
  #This function works by constructing a marssMODEL object (form=marss) with parameter values at each t to be predicted
  #This will be passed to the Kalman filter to compute the expected x values
  #Then to haty to get the expected y values
  
  #newdata must reflect the inputs allowed for a specific form; y (data) is an input.
  #n.ahead is the number of time steps
  #t.start is where to start relative to the original data.  
  #t.start=1 would start at the beginning of the original data
  #the initial x will be E(x)_(t.start-1)
  
  if(class(object)!="marssMLE") stop("predict.marssMLE: object must be a marssMLE object.",call.=FALSE)
  
  #if n.ahead is not passed in AND newdata is passed in, n.ahead will be taken from the info in newdata
  if(missing(n.ahead) & length(newdata)!=0) n.ahead=NULL
  #By default, start the predict at the next time step after the end of the time series
  if(missing(t.start)) t.start=dim(object$marss$data)[2]+1
  #make sure t.start is not bigger than TT; that doesn't make sense
  if(!all(is.wholenumber(t.start))) stop("predict.marssMLE: t.start must be whole numbers.",call.=FALSE)  
  if(max(t.start)>(dim(object$marss$data)[2]+1)) stop("predict.marssMLE: t.start cannot be greater than length of data+1.  Leave off to start where the data ends.",call.=FALSE)
  if(min(t.start)<1) stop("predict.marssMLE: t.start must be between 1 and the length of data.  Leave off to start at the first time step after the data ends.",call.=FALSE)
  
  rtn.list=list()
  orig.obj=object
  
  #By default newdata is assumed to be in a form compatible with $call$form
  if(is.null(object[["call"]][["form"]])){ 
    form=attr(object[["model"]][["form"]],"form") }else{ form=object[["call"]][["form"]] }
  
  #First make sure specified equation form has a corresponding function to do the conversion to marssMODEL object
  predict.fun = paste("predict_",form[1],sep="")
  tmp=try(exists(predict.fun,mode="function"),silent=TRUE)
  if(isTRUE(tmp)){
    #the predict function takes the marssMLE obj along with newdata and returns a
    #marssMODEL (form=marss) object constructed using info in newdata that is ready for use in prediction
    #x means the x (marssMLE object) for prediction
    #also returns the newdata list
    tmp=eval(call(predict.fun, object, newdata, n.ahead, t.start))
    pred.obj=tmp$MLEobj
    newdata=tmp$newdata
  }else{ stop(paste("predict.marssMLE: missing",predict.fun,"function to interpret newdata."), call.=FALSE) }
  
  marss.model=pred.obj[["marss"]]
  ## Check that the marssMODEL object is ok
  ## More checking on the control list is done by is.marssMLE() to make sure the MLEobj is ready for fitting
  tmp = is.marssMODEL(marss.model, method=object[["method"]])
  
  if(!isTRUE(tmp)) {
    cat("predict_marssMLE: Stopped in predict() due to problem(s) with model specification. \nThe marssMLE object constructed for predict call is being returned (invisibly).  Examine $model for problems.\n")
    cat(tmp)
    return(invisible(pred.obj))
  }
  kf=MARSSkf(pred.obj)
  
  Ey=MARSShatyt(pred.obj)  #this needs the $kf element set and $data
  rtn.list$E.x=kf$xtT
  if(se.fit & !is.null(kf[["VtT"]])){
    m = dim(pred.obj$marss$fixed$x0)[1]
    TT = dim(pred.obj$marss$data)[2]
    if(m == 1) states.se = sqrt(matrix(kf$VtT[,,1:TT], nrow=1))
    if(m > 1) {
      states.se = matrix(0, nrow=m, ncol=TT)
      for(i in 1:TT) 
        states.se[,i] = t(sqrt(takediag(kf$VtT[,,i])))
    }
  }else{  states.se=NULL }
  rtn.list$x.se = states.se
  
  rtn.list$E.y=Ey[["ytT"]]
  if(se.fit & !is.null(Ey[["OtT"]])){
    n = dim(pred.obj$marss$data)[1]
    TT = dim(pred.obj$marss$data)[2]
    var.y.tT=
      if(n == 1){
        var.y.tT = Ey[["OtT"]][,,1:TT] - Ey[["ytT"]][,1:TT]*Ey[["ytT"]][,1:TT]
        y.se = sqrt(matrix(var.y.tT, nrow=1))
      }
    if(n > 1) {
      y.se = matrix(0, nrow=n, ncol=TT)
      for(i in 1:TT){
        var.y.tT = Ey[["OtT"]][,,i] - Ey[["ytT"]][,i,drop=FALSE]%*%t(Ey[["ytT"]][,i,drop=FALSE])
        y.se[,i] = t(sqrt(takediag( var.y.tT )))
      }
    }
  }else{  y.se=NULL }
  rtn.list$y.se = y.se
  
  if(prediction.intervals){
    pred.obj$kf=NULL #just to be sure
    pred.obj$Ey=NULL #just to be sure
    #Calculation of prediction intervals
    
    #generate a bootstrapped parameter set from original model
    par.boot=MARSSboot(orig.obj, nboot=nboot, output="parameters", param.gen=param.gen, control=NULL, silent=FALSE)$boot.params
    #need to get rid of any x0 and V0 related parameters since these are not used
    tmp.obj=orig.obj; tmp.obj$par$x0=matrix(0,0,1); tmp.obj$par$V0=matrix(0,0,1)
    ok=rownames(par.boot) %in% names(MARSSvectorizeparam(tmp.obj))
    states.boot=array(NA,dim=c(attr(pred.obj$marss,"model.dims")[["x"]],nboot))
    y.boot=array(NA,dim=c(attr(pred.obj$marss,"model.dims")[["y"]],nboot))
    model.dims=attr(pred.obj$marss,"model.dims")
    
    for(i in 1:nboot){
      #Take those par boot and put in the MLEobj for prediction (the one made with newdata)
      MLEobj.new.boot=MARSSvectorizeparam(pred.obj,par.boot[ok,i])
      
      #compute the initial xtT at t.start-1 for these values
      tmp.obj=MARSSvectorizeparam(orig.obj,par.boot[,i])
      #Set the initial conditons
      #Set these to the estimated distribution of x at t.start-1 conditioned on all the data (original)
      kf=MARSSkf(tmp.obj)
      if(t.start==1){
        MLEobj.new.boot$marss$fixed$x0=array(kf$x0T, dim=model.dims$x0)
        MLEobj.new.boot$marss$fixed$V0=array(vec(kf$V0T), dim=c(model.dims$V0[1]*model.dims$V0[2],1,1))
      }else{
        MLEobj.new.boot$marss$fixed$x0=array(kf$xtT[,t.start-1], dim=model.dims$x0)
        MLEobj.new.boot$marss$fixed$V0=array(vec(kf$VtT[,,t=(t.start-1)]), dim=c(model.dims$V0[1]*model.dims$V0[2],1,1))
      }
      
      states.boot[,,i]=MARSSkf(MLEobj.new.boot)$xtT
      y.boot[,,i]=MARSShatyt(MLEobj.new.boot)$ytT
    }
  }
  if(verbose){
    rtn.list$states.boot=states.boot
    rtn.list$y.boot=y.boot
    rtn.list$par.boot=par.boot
  }
  
  rtn.list$newdata=newdata
  rtn.list$pred.MLEobj=pred.obj
  invisible(rtn.list)
}  #end of predict.marssMLE
