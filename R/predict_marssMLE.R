######################################################################################################  predict method for class marssMLE. Prediction intervals
##################################################################################
predict.marssMLE <- function(object, h=0,
                              conf.level = c(0.80, 0.95),
                              type = c("ytT", "xtT"),
                              newdata = list(t=NULL, y=NULL, c=NULL, d=NULL),
                              interval = c("prediction", "confidence", "none"),
                              fun.kf = c("MARSSkfas", "MARSSkfss"),
                              x0 = list(x0=NULL, tinitx=NULL),
                              ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  fun.kf <- match.arg(fun.kf)
  if(object[["fun.kf"]] != fun.kf) message(paste0(fun.kf, "is being used for prediction. This is different than fun.kf in the marssMLE object.\n"))
  MODELobj <- object[["model"]]
  model.dims <- attr(MODELobj, "model.dims")
  TT <- model.dims[["y"]][2]
  nx <- switch(type,
               ytT = model.dims[["y"]][1],
               xtT = model.dims[["x"]][1])
  
  if (is.null(object[["par"]])) {
    stop("predict.marssMLE: The marssMLE object does not have the par element.  Most likely the model has not been fit.\n", call. = FALSE)
  }
  if( interval == "none" ) conf.level <- c()
  if (length(conf.level) > 0 && (!is.numeric(conf.level) ||  any(conf.level > 1) || any(conf.level < 0)))
    stop("predict.marssMLE: conf.level must be between 0 and 1.", call. = FALSE)
  if(length(conf.level) == 0) interval <- "none"
  if (!missing(h) && (!is.numeric(h) ||  length(h) != 1 || h < 0 || (h %% 1) != 0))
    stop("predict.marssMLE: h must be an integer > 0.\n", call. = FALSE)
  if (!is.list(x0))
    stop("predict.marssMLE: x0 must be list with elements x0 and/or tinitx.\n", call. = FALSE)
  if (!is.null(x0$tinitx) && is.null(x0$x0))
    stop("predict.marssMLE: if x0$tinitx is specified, so must x0$x0.\n", call. = FALSE)
  if(missing(x0) || all(unlist(lapply(x0, is.null)))) x0.estimate <- FALSE else x0.estimate <- TRUE
  extras <- list()
  if (!missing(...)) {
    extras <- list(...)
  }
  
  if(!missing(h)){
    if(!missing(x0)) message("predict.marssMLE: x0 was passed in. Will be ignored since h also passed in.\n")

    outlist <- forecast.marssMLE(object, h=h, conf.level = conf.level,
                                       interval = interval,
                                       type = type, newdata = newdata,
                                       fun.kf = fun.kf, ...)
    tmp <- colnames(outlist[["pred"]])
    colnames(outlist[["pred"]])[which(tmp==".sd")] <- "se"
    outlist$x0 <- coef(object, type="matrix")[["x0"]]
    outlist$tinitx <- object[["model"]][["tinitx"]]
    
    class(outlist) <- "marssPredict"
    
    return(outlist)
  }
  
  # h not passed in, start checking newdata
  
  if(!is.null(newdata[["t"]])){
    if(!is.numeric(newdata[["t"]]) || !is.vector(newdata[["t"]]))
      stop("predict.marssMLE: t in newdata must be a numeric vector.\n", call.=FALSE)
    if(!all(diff(newdata[["t"]])==1) || newdata[["t"]][1] < 1)
      stop("predict.marssMLE: t in newdata must be a positive ordered sequence of integers, one time step apart (like 1,2,3...).\n", call.=FALSE)
  }
  nonewdata <- all(unlist(lapply(newdata[c("y", "c", "d")], is.null)))
  if(nonewdata && !is.null(newdata[["t"]])){
    if(length(newdata[["t"]]) > TT || newdata[["t"]][length(newdata[["t"]])] > TT)
      stop("predict.marssMLE: if no y, c or d in newdata, t must be within the original time steps.\n", call.=FALSE)
  }
  if(nonewdata && is.null(newdata[["t"]])) newdata[["t"]] <- 1:TT
  
  # h is NULL if here; if newdata is missing and so is h, then use model data
  if(nonewdata){
    if(!missing(x0)) message("predict.marssMLE: x0 was passed in. Will be ignored since newdata was not passed in.\n")
    
    newMLEobj <- object
    
  } else {
  
    # newdata was passed in. Need to make newMLEobj
  if(!is.list(newdata) || !all(names(newdata) %in% c("t", "y", "c", "d")))
    stop("predict.marssMLE: newdata must be list with only t, y, c and/or d.\n", call.=FALSE)
  
  # We need the model in marxss form
  new.MODELlist <- coef(object, type="matrix", form="marxss")
  
  isxreg <- list(c=TRUE, d=TRUE, y=TRUE)
  for(elem in c("c", "d")){
    tmp <- coef(object, type="matrix", form="marxss")[[elem]]
    if(dim(tmp)[1]==1 && dim(tmp)[2]==1 && tmp==0){
      isxreg[[elem]] <- FALSE
    }else{
      dim(new.MODELlist[[elem]]) <- model.dims[[elem]][c(1,3)]
    }
  }
  for(elem in c("y", "c", "d")){
    if(!isxreg[[elem]]){
      if(!is.null(newdata[[elem]])) message(paste0("predict.marssMLE(): model does not include ", elem, ". ", elem, " in newdata is being ignored.\n"))
    }else{
      if(is.null(newdata[[elem]])){ 
        stop(paste0("predict.marssMLE(): model includes ", elem, ". ", elem, " must be in newdata.\n"), call.=FALSE)
      }
      if(!is.numeric(newdata[[elem]])){ 
        stop("predict.marssMLE: y, c, and d in newdata must be numeric (use class() and is.numeric() to test what you are passing in).\n", call.=FALSE)
      }
      if(is.vector(newdata[[elem]])) newdata[[elem]] <- matrix(newdata[[elem]], nrow = 1)
      if(!is.matrix(newdata[[elem]])) stop(paste0("predict.marssMLE(): newdata ", elem, " must be a matrix with ", model.dims[[elem]][1], " rows.\n"), call.=FALSE)
      if(dim(newdata[[elem]])[1]!=model.dims[[elem]][1]) stop(paste0("predict.marssMLE(): model ", elem, " has ", model.dims[[elem]][1], " rows. ", elem, " in newdata does not.\n"), call.=FALSE)

      if(elem != "y") new.MODELlist[[elem]] <- newdata[[elem]]
      if(elem == "y") new.data <- newdata[[elem]]
    }
  }
  # check that ncols match
  tmp <- unlist(lapply(newdata[c("y", "c", "d")], ncol))
  ncol.newdata <- tmp[1]
  if(!all(tmp==ncol.newdata))  stop("predict.marssMLE(): y, c, and d in newdata must all have the same number of columns.\n", call.=FALSE)
  # check that t length matches ncols
  if(!is.null(newdata[["t"]]) && length(newdata[["t"]])!=ncol.newdata)
    stop("predict.marssMLE(): t in newdata must be the same length as the number of columns in y, c and d.\n", call.=FALSE)
  if(!x0.estimate)
    message("predict.marssMLE(): x0 from the original data is being used for prediction.\n")
  
  # check if parameters are time-varying. If so, t used to specify which parameters to use
  for(elem in names(new.MODELlist)){
    if(elem %in% c("c", "d", "x0", "V0")) next
    if(model.dims[[elem]][3] == TT){
      if(newdata[["t"]][length(newdata[["t"]])] > TT) stop("predict.marssMLE(): the model has time-varying parameters. t in newdata needed in order to specify which parameters to use.\n", call.=FALSE)
      tmp <- array(NA, dim=c(model.dims[[elem]][1:2], model.dims[[elem]][3]+h))
      tmp[,,1:TT] <- new.MODELlist[[elem]]
      tmp[,,(TT+1):(TT+h)] <- new.MODELlist[[elem]][,,TT, drop=FALSE]
      new.MODELlist[[elem]] <- tmp
      message(paste0(elem, " is time-varying. The value at ", TT, " is used for the forecast.\n"))
    }
  }
  #Passed all checks. Can set t now
  if(is.null(newdata[["t"]])) newdata[["t"]] <- 1:ncol.newdata

  
  if(x0.estimate){
    new.MODELlist[["tinitx"]] <- x0[["tinitx"]]
    new.MODELlist[["x0"]] <- x0[["x0"]]
    newMLEobj <- MARSS(newdata[["y"]], model=new.MODELlist, silent=TRUE, method=object[["method"]])
  }
  
  if(!x0.estimate){
  new.MODELlist[["tinitx"]] <- object[["model"]][["tinitx"]]
  newMLEobj <- MARSS(newdata[["y"]], model=new.MODELlist, silent=TRUE, method=object[["method"]])
  }
} # end setting up newMLEobj
  
  if(type=="ytT"){
    if(length(conf.level) == 0) interval <- "none"
    cols <- switch(interval,
                   prediction = c(".rownames", "t", "y", ".fitted", ".sd.y", ".lwr", ".upr"),
                   none = c(".rownames", "t", "y", ".fitted"),
                   confidence = c(".rownames", "t", "y", ".fitted", ".se.fit", ".conf.low", ".conf.up"))
    ret <- fitted.marssMLE(newMLEobj, type=type, interval=interval, conf.level=conf.level[1])[cols]
    colnames(ret)[which(colnames(ret)==".fitted")] <- "estimate"
    colnames(ret)[which(colnames(ret)==".sd.y")] <- "se"
    colnames(ret)[which(colnames(ret)==".se.fit")] <- "se"
    colnames(ret)[which(colnames(ret)==".lwr")] <- paste("Lo", 100*conf.level[1])
    colnames(ret)[which(colnames(ret)==".upr")] <- paste("Hi", 100*conf.level[1])
    colnames(ret)[which(colnames(ret)==".conf.low")] <- paste("Lo", 100*conf.level[1])
    colnames(ret)[which(colnames(ret)==".conf.up")] <- paste("Hi", 100*conf.level[1])
    if(interval != "none" && length(conf.level) > 1){
      for(i in 2:length(conf.level)){
        cols <- switch(interval,
                       prediction = c(".lwr", ".upr"),
                       confidence = c(".conf.low", ".conf.up"))
        tmp <- fitted.marssMLE(newMLEobj, type=type, interval=interval, conf.level=conf.level[i])[cols]
        colnames(tmp) <- paste(c("Lo", "Hi"), 100*conf.level[i])
        ret <- cbind(ret, tmp)
      }
    }
  }
  if(type=="xtT"){
    conf.int <- TRUE
    if(length(conf.level) == 0 || interval == "none") conf.int <- FALSE
    ret <- tidy.marssMLE(newMLEobj, type=type, conf.int=conf.int, conf.level=conf.level[1])
    colnames(ret)[which(colnames(ret)=="std.error")] <- "se"
    colnames(ret)[which(colnames(ret)=="conf.low")] <- paste("Lo", 100*conf.level[1])
    colnames(ret)[which(colnames(ret)=="conf.high")] <- paste("Hi", 100*conf.level[1])
    if(conf.int && length(conf.level) > 1){
      for(i in 2:length(conf.level)){
        tmp <- tidy.marssMLE(newMLEobj, type=type, conf.int=conf.int, conf.level=conf.level[i])[c("conf.low", "conf.high")]
        colnames(tmp) <- paste(c("Lo", "Hi"), 100*conf.level[i])
        ret <- cbind(ret, tmp)
      }
    }
  }
  
  if(nonewdata){
    ret <- subset(ret, t %in% newdata[["t"]])
  }
  
  # set t in ret with t in newdata
  ret$t <- rep(newdata[["t"]], nx)
  
  # Set up output
  outlist <- list(
    method = c("MARSS", object[["method"]]),
    model = object,
    newdata = newdata,
    level = 100*conf.level,
    pred = ret,
    type = type,
    t = newdata[["t"]],
    h = h,
    x0 = coef(newMLEobj, type="matrix")[["x0"]],
    tinitx = newMLEobj[["model"]][["tinitx"]]
  )
  
  class(outlist) <- "marssPredict"
  
  return(outlist)
  
} 
