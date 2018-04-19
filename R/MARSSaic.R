#######################################################################################################
#   MARSSaic function
#   Adds to marssMLE object:
#   elements in output arg
#   samp.size, num.params
#######################################################################################################
MARSSaic = function(MLEobj, output=c("AIC","AICc"), Options=list(nboot=1000, return.logL.star=FALSE,
  silent=FALSE)) {
  # Options$nboot is the number of bootstrap replicates to do
  # mssm.model is a specified model (needs the structure and parameter elements of the list)
  # output tells what output to produce; this can be a vector  if multiple items should be returned
  #      AIC, AICc, AICbp, AICbb, AICi, boot.params
  # silent is a flag to indicate whether the progress bar should be printed

  if(!is.list(Options)){
    msg = " Options argument must be passed in as a list.\n"
    cat("\nErrors were caught in MARSSaic \n", msg, sep="") 
    stop("Stopped in MARSSaic() due to problem(s) with arguments.\n", call.=FALSE)
  } 
  ## Set options if some not passed in
  if(is.null(Options[["nboot"]])) Options$nboot = 1000
  if(is.null(Options[["return.logL.star"]])) Options$return.logL.star = FALSE
  if(is.null(Options[["silent"]])) Options$silent = FALSE
  if(class(MLEobj)[1]!="marssMLE") {
    stop("Stopped in MARSSaic(). An object of class marssMLE is required.\n", call.=FALSE)
    }
    
  if(is.null(MLEobj[["logLik"]])){
    msg = " No log likelihood.  This function expects a model fitted via maximum-likelihood.\n"
    cat("\n","Errors were caught in MARSSaic \n", msg, sep="") 
    stop("Stopped in MARSSaic() due to problem(s) with the MLE object passed in.\n", call.=FALSE)
    }
  return.list=list()

  ## Some renaming for readability
  model = MLEobj$marss
  #kf = MLEobj$kf
  loglike = MLEobj$logLik
  
  ##### AIC and AICc calculations
  if("AIC" %in% output | "AICc" %in% output){
    K = 0
    for (elem in c("Z", "A", "B", "U", "x0", "R", "Q","V0" )) {
      pars = dim(model$free[[elem]])[2]
      K = K + pars
    }
    MLEobj$AIC = -2*loglike + 2*K
    samp.size = sum(!is.na(model$data))
    MLEobj$AICc = ifelse(samp.size > (K+1),-2*loglike + 2*K*(samp.size/(samp.size-K-1)),"NA, number of data points less than K+1")
    MLEobj$samp.size = samp.size 
    MLEobj$num.params = K
  }

  ##### AICbb & AICbp  
  method = c("AICbb","AICbp") 

  for (m in method[which(method %in% output)]) {
    bootstrap.method = switch(m,
			AICbb = "innovations",
			AICbp = "parametric")
    logL.star=0

    drawProgressBar = FALSE #If the time library is not installed, no prog bar
    if(!Options$silent) { #then we can draw a progress bar
       cat(paste(m, "calculation in progress...\n"))
       prev = progressBar()
       drawProgressBar = TRUE
    }

    boot.params = MARSSboot(MLEobj, nboot=Options$nboot, output="parameters", sim=bootstrap.method,
          param.gen="MLE", silent=TRUE)$boot.params

    if((bootstrap.method == "parametric") && ("boot.params" %in% output)) MLEobj$boot.params = boot.params

    for(i in 1:Options$nboot) {
      boot.model = MARSSvectorizeparam( MLEobj, parvec=boot.params[,i] )  #boot.model is a MLEobj
      logL.star[i] = MARSSkf(boot.model, only.logLik=TRUE)$logLik	
      if(drawProgressBar) prev = progressBar(i/Options$nboot,prev)
    }
    MLEobj[[m]] = -4*(sum(logL.star))/Options$nboot + 2*loglike   # -2*model$loglik + 2*(1/N)*(-2)*sum(boot$Yloglike-model$logLik)
    if(Options$return.logL.star == TRUE) {
      tmp = paste(m, ".logL.star", sep="")
      MLEobj[[tmp]] = logL.star
    }
  }

return(MLEobj)
}
