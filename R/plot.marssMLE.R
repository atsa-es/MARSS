plot.marssMLE <- 
  function(x, 
           plot.type=c("observations", "states", "model.residuals", "state.residuals", "model.residuals.qqplot", "state.residuals.qqplot"), 
           form=c("marxss", "marss", "dfa"), 
           conf.int=TRUE, conf.level=0.95, decorate=TRUE, ...)
  {
    
    # Argument checks
    plot.type = match.arg(plot.type, several.ok = TRUE)
    if(!is.numeric(conf.level) | conf.level>1 | conf.level < 0) stop("plot.marssMLE: conf.level must be between 0 and 1.", call. = FALSE)
    if(!(conf.int%in%c(TRUE,FALSE))) stop("plot.marssMLE: conf.int must be TRUE/FALSE", call. = FALSE)
    
    if(missing(form)){
      model_form = attr(x[["model"]],"form")[1]
    }else{
      model_form = match.arg(form)
    }
    
    extras = list()
    if(!missing(...)){
      extras=list(...)
      allowednames=c("rotate", "method", "hessian.fun", "nboot")
      bad.names=names(extras)[!(names(extras)%in%allowednames)]
      if(!all(names(extras)%in%allowednames)) stop(paste("plot.marssMLE:", paste(bad.names, collapse=" "), "is/are unknown argument(s). See ?tidy.marssMLE for allowed arguments.\n"), call. = FALSE)
    }
    if(model_form!="dfa" & "rotate"%in%names(extras)){
      cat("plot.marssMLE: 'rotate' argument is ignored if form!='dfa'\n Pass in form='dfa' if your model is a DFA model, but the form \n attribute is not set (because you set up your DFA model manually).\n\n")
      rotate=FALSE
    }
    # End Argument checks
    
    alpha = 1-conf.level
    
    if("states" %in% plot.type) {
      # make plot of states and CIs
      
      if("rotate"%in%names(extras)){
        rotate=extras[["rotate"]]
        if(!(rotate %in% c(TRUE, FALSE))) stop("tidy.marssMLE: rotate must be TRUE/FALSE. \n")
      }else{ rotate=FALSE }
      
      states = tidy.marssMLE(x, type="states", conf.int=conf.int, conf.level=conf.level, ...)
      if(model_form == "dfa"){
        if(rotate){ rottext = "rotated" }else{ rottext = "" }
        states$term = paste("DFA", rottext, "trend",states$term)
      }else{
        states$term = paste0("State ",states$term)
      }
      
      nX = min(9,attr(x$model, "model.dims")$x[1])
      plot.nrow = round(sqrt(nX))
      plot.ncol = ceiling(nX/plot.nrow)
      par(mfrow=c(plot.nrow, plot.ncol), mar=c(2, 4, 2, 1) + 0.1)
      for(plt in unique(states$term)){
        with(subset(states, states$term==plt), {
          ylims=c(min(estimate,conf.low,na.rm=TRUE),max(estimate,conf.high,na.rm=TRUE))
          plot(t, estimate, type="l", xlab="", ylab="Estimate", ylim=ylims)
          title(plt)
          if(conf.int) polygon(c(t, rev(t)), c(conf.low, rev(conf.high)),col="grey",border=FALSE)
          lines(t,estimate)
          box()
        })
      }
      plot.type = plot.type[plot.type!="states"]
      cat(paste("plot type = Estimated States\n"))
      if(length(plot.type) != 0){
        ans <- readline(prompt="Hit <Return> to see next plot (q to exit): ")
        if(tolower(ans)=="q") return()
      }
    }
    
    if("observations" %in% plot.type) {
      # make plot of observations
      df = augment.marssMLE(x, "observations", form=model_form)
      df$ymin = df$.fitted - qnorm(alpha/2)*df$.se.fit
      df$ymax = df$.fitted + qnorm(alpha/2)*df$.se.fit
      nY = min(9,attr(x$model, "model.dims")$y[1])
      plot.ncol = round(sqrt(nY))
      plot.nrow = ceiling(nY/plot.ncol)
      par(mfrow=c(plot.nrow, plot.ncol), mar=c(2, 4, 2, 1) + 0.1)
      for(plt in levels(df$.rownames)){
        with(subset(df, .rownames==plt), {
          ylims=c(min(.fitted,y,ymin,ymax,na.rm=TRUE),max(.fitted,y,ymin,ymax,na.rm=TRUE))
          plot(t,.fitted, type="l",xlab="",ylab="Estimate",ylim=ylims)
          title(plt)
          if(conf.int) polygon(c(t, rev(t)), c(ymin, rev(ymax)),col="grey",border=FALSE)
          lines(t,.fitted)
          points(t,y,col="blue",pch=19)
          box()
        })
      }
      plot.type = plot.type[plot.type!="observations"]
      cat(paste("plot type = Observations with Fitted Values\n"))
      if(length(plot.type) != 0){
        ans <- readline(prompt="Hit <Return> to see next plot (q to exit): ")
        if(tolower(ans)=="q") return()
      }
    }
    
    if("model.residuals" %in% plot.type) {
      # make plot of observation residuals
      df = augment.marssMLE(x, "observations", form="marxss")
      df$.resids[is.na(df$y)]=NA
      nY = min(9,attr(x$model, "model.dims")$y[1])
      plot.ncol = round(sqrt(nY))
      plot.nrow = ceiling(nY/plot.ncol)
      par(mfrow=c(plot.nrow, plot.ncol), mar=c(2, 4, 2, 1) + 0.1)
      for(plt in levels(df$.rownames)){
        with(subset(df, .rownames==plt), {
          ylims=c(min(.resids,na.rm=TRUE),max(.resids,na.rm=TRUE))
          if(decorate){ 
            lo = predict(loess(.resids~t), newdata=data.frame(t=t),se=TRUE)
            lo.t=names(lo$fit)
            alp=qnorm((1-conf.level)/2)*lo$se.fit
            ymin=alp + lo$fit
            ymax=-1*alp + lo$fit
            ylims=c(min(.resids,ymin,na.rm=TRUE),max(.resids,ymax,na.rm=TRUE))
          }
          plot(t,.resids, type="p", xlab="",
               ylab="",col="blue",ylim=ylims)
          title(plt)
          if(decorate){ 
            polygon(c(t, rev(t)), 
                    c(ymin, rev(ymax)),
                    col="grey",border=FALSE)
            lines(t,lo$fit, col="blue")
          }
          points(t,.resids,col="blue",pch=19)
          box()
          abline(h=0, lty=3)
        })
        mtext("Observation residuals, y - E[y]",side=2,outer=TRUE,line=-1)
      }
      plot.type = plot.type[plot.type!="model.residuals"]
      cat(paste("plot type = Model Residuals\n"))
      if(length(plot.type) != 0){
        ans <- readline(prompt="Hit <Return> to see next plot (q to exit): ")
        if(tolower(ans)=="q") return()
      }
    }
    
    if("state.residuals" %in% plot.type) {
      # make plot of process residuals; set form='marxss' to get process resids
      df = augment.marssMLE(x, "states", form="marxss")
      df$.rownames = paste0("State ",df$.rownames)
      nX = min(9,attr(x$model, "model.dims")$x[1])
      plot.nrow = round(sqrt(nX))
      plot.ncol = ceiling(nX/plot.nrow)
      par(mfrow=c(plot.nrow, plot.ncol), mar=c(2, 4, 2, 1) + 0.1)
      for(plt in unique(df$.rownames)){
        with(subset(df, .rownames==plt), {
          ylims=c(min(.resids,na.rm=TRUE),max(.resids,na.rm=TRUE))
          if(decorate){ 
            lo = predict(loess(.resids~t), newdata=data.frame(t=t),se=TRUE)
            lo.t=names(lo$fit)
            alp=qnorm((1-conf.level)/2)*lo$se.fit
            ymin=alp + lo$fit
            ymax=-1*alp + lo$fit
            ylims=c(min(.resids,ymin,na.rm=TRUE),max(.resids,ymax,na.rm=TRUE))
          }
          plot(t,.resids, type="p", xlab="",
               ylab="",col="blue",ylim=ylims)
          title(plt)
          if(decorate){ 
            polygon(c(t, rev(t)), 
                    c(ymin, rev(ymax)),
                    col="grey",border=FALSE)
            lines(t,lo$fit, col="blue")
          }
          points(t,.resids,col="blue",pch=19)
          box()
          abline(h=0, lty=3)
        })
        mtext("State residuals, xtT - E[x]",side=2,outer=TRUE,line=-1)
      }
      plot.type = plot.type[plot.type!="state.residuals"]
      cat(paste("plot type = State Residuals\n"))
      if(length(plot.type) != 0){
        ans <- readline(prompt="Hit <Return> to see next plot (q to exit): ")
        if(tolower(ans)=="q") return()
      }
    }
    
    slp = function(yy){
      y <- quantile(yy[!is.na(yy)], c(0.25, 0.75))
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y)/diff(x)
      return(slope)
    }
    int = function(yy){
      y <- quantile(yy[!is.na(yy)], c(0.25, 0.75))
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y)/diff(x)
      int <- y[1L] - slope * x[1L]
      return(int)
    }
    
    if("model.residuals.qqplot" %in% plot.type) {
      # make plot of observation residuals
      df = augment.marssMLE(x, "observations", form="marxss")
      slope=tapply(df$.std.resid,df$.rownames,slp)
      intercept=tapply(df$.std.resid,df$.rownames,int)
      nY = min(9,attr(x$model, "model.dims")$y[1])
      plot.ncol = round(sqrt(nY))
      plot.nrow = ceiling(nY/plot.ncol)
      par(mfrow=c(plot.nrow, plot.ncol), mar=c(2, 4, 2, 1) + 0.1)
      for(plt in levels(df$.rownames)){
        with(subset(df, df$.rownames==plt), {
          qqnorm(.std.resid, main=plt)
          abline(a=intercept[plt], b=slope[plt], col="blue")
        })
      }
      plot.type = plot.type[plot.type!="model.residuals.qqplot"]
      cat(paste("plot type = Model Standardized Residuals\n"))
      if(length(plot.type) != 0){
        ans <- readline(prompt="Hit <Return> to see next plot (q to exit): ")
        if(tolower(ans)=="q") return()
      }
    }  
    
    if("state.residuals.qqplot" %in% plot.type) {
      # make qqplot of state residuals
      df = augment.marssMLE(x, "states", form="marxss")
      df$.rownames = paste0("State ",df$.rownames)
      slope=tapply(df$.std.resid,df$.rownames,slp)
      intercept=tapply(df$.std.resid,df$.rownames,int)
      nX = min(9,attr(x$model, "model.dims")$x[1])
      plot.nrow = round(sqrt(nX))
      plot.ncol = ceiling(nX/plot.nrow)
      par(mfrow=c(plot.nrow, plot.ncol), mar=c(2, 4, 2, 1) + 0.1)
      for(plt in unique(df$.rownames)){
        with(subset(df, .rownames==plt), {
          qqnorm(.std.resid, main=plt)
          abline(a=intercept[plt], b=slope[plt], col="blue")
        })
      }
      plot.type = plot.type[plot.type!="state.residuals.qqplot"]
      cat(paste("plot type = Standardized State Residuals\n"))
      if(length(plot.type) != 0){
        ans <- readline(prompt="Hit <Return> to see next plot (q to exit): ")
        if(tolower(ans)=="q") return()
      }
    }
    
  }