plot.marssMLE <- 
  function(x, 
           plot.type=c("observations", "states", "model.residuals", "state.residuals", "model.residuals.qqnorm", "state.residuals.qqnorm"), 
           form=c("marxss", "marss", "dfa"), ...) {
  plot.type = match.arg(plot.type, several.ok = TRUE)
  if(missing(form)){
    model_form = attr(x[["model"]],"form")
  }else{
    form = match.arg(form)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for plot.marrsMLE. Please install it.", call. = FALSE)
  }  
  plts = list()
  
  if("states" %in% plot.type) {
    # make plot of states and CIs
    states = tidy.marssMLE(x, "states")
    if(model_form[1] == "dfa"){
      states$term = paste0("DFA trend ",states$term)
    }else{
      states$term = paste0("State ",states$term)
    }
    p1 = ggplot2::ggplot(data=states, ggplot2::aes_(~t, ~estimate)) +
      ggplot2::geom_ribbon(data=states, ggplot2::aes_(ymin = ~conf.low, ymax = ~conf.high), alpha=0.3, col="grey") +
      ggplot2::geom_line() + 
      ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
      ggplot2::facet_wrap(~term, scale="free_y")
    plts[["states"]] = p1
    if(identical(plot.type, "states")) return(p1)
  }
  
  if("observations" %in% plot.type) {
    # make plot of observations
    df = augment.marssMLE(x, "observations")
    df$ymin = df$.fitted - 1.96*df$.se.fit
    df$ymax = df$.fitted + 1.96*df$.se.fit
    p2 = ggplot2::ggplot(data=df, ggplot2::aes_(~t, ~.fitted)) +
      ggplot2::geom_ribbon(data=df, ggplot2::aes_(ymin = ~ymin, ymax= ~ymax), alpha=0.3, col="grey") +
      ggplot2::geom_line() + 
      ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
      ggplot2::facet_wrap(~.rownames, scale="free_y") + 
      ggplot2::geom_point(data=df[!is.na(df$y),], ggplot2::aes_(~t, ~y), col="blue")
    plts[["observations"]] = p2
    if(identical(plot.type, "observations")) return(p2)
  }
  
  if("model.residuals" %in% plot.type) {
    # make plot of observation residuals
    df = augment.marssMLE(x, "observations")
    p1 = ggplot2::ggplot(df[(!is.na(df$.resids) & !is.na(df$y)),], ggplot2::aes_(~t, ~.resids)) + 
      ggplot2::geom_point(col="blue") +
      ggplot2::stat_smooth(method="loess") +
      ggplot2::xlab("Time") + 
      ggplot2::ylab("Observation residuals, y - E[y]") +
      ggplot2::facet_wrap(~.rownames, scale="free_y") + 
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype=3)
    plts[["model.residuals"]] = p1
    if(identical(plot.type, "model.residuals")) return(p1)
  }  
  
  if("state.residuals" %in% plot.type) {
    # make plot of process residuals
    df = augment.marssMLE(x, "states")
    df$.rownames = paste0("State ",df$.rownames)
    p1 = ggplot2::ggplot(df[!is.na(df$.resids),], ggplot2::aes_(~t, ~.resids)) + 
      ggplot2::geom_point(col="blue") + 
      ggplot2::stat_smooth(method="loess") +
      ggplot2::xlab("Time") + 
      ggplot2::ylab("State residuals, xtT - E[x]") +
      ggplot2::facet_wrap(~.rownames, scale="free_y") + 
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype=3)
    plts[["state.residuals"]] = p1
    if(identical(plot.type, "state.residuals")) return(p1)
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
  
  if("model.residuals.qqnorm" %in% plot.type) {
    # make plot of observation residuals
    df = augment.marssMLE(x, "observations")
    slope=tapply(df$.std.resid,df$.rownames,slp)
    intercept=tapply(df$.std.resid,df$.rownames,int)
    abline.dat=data.frame(.rownames=names(slope), slope=slope, intercept=intercept)
    p1 = ggplot2::ggplot(df) +
      geom_qq(aes_(sample = ~.std.resid),na.rm=TRUE) +
      ggplot2::xlab("Theoretical Quantiles") + 
      ggplot2::ylab("Standardized Model Residuals") +
      ggplot2::geom_abline(data=abline.dat, ggplot2::aes_(slope=~slope, intercept=~intercept),color="blue") +
      ggplot2::facet_wrap(~.rownames, scale="free_y")
    plts[["model.residuals.qqnorm"]] = p1
    if(identical(plot.type, "model.residuals.qqnorm")) return(p1)
  }  
  
  if("state.residuals.qqnorm" %in% plot.type) {
    # make qqplot of state residuals
    df = augment.marssMLE(x, "states")
    df$.rownames = paste0("State ",df$.rownames)
    slope=tapply(df$.std.resid,df$.rownames,slp)
    intercept=tapply(df$.std.resid,df$.rownames,int)
    abline.dat=data.frame(.rownames=names(slope), slope=slope, intercept=intercept)
    p1 = ggplot2::ggplot(df) +
      geom_qq(aes_(sample = ~.std.resid),na.rm=TRUE) +
      ggplot2::xlab("Theoretical Quantiles") + 
      ggplot2::ylab("Standardized State Residuals") +
      ggplot2::geom_abline(data=abline.dat, ggplot2::aes_(slope=~slope, intercept=~intercept),color="blue") +
      ggplot2::facet_wrap(~.rownames, scales="free_y")
    plts[["state.residuals.qqnorm"]] = p1
    if(identical(plot.type, "state.residuals.qqnorm")) return(p1)
  }    
  for(i in plot.type){
    print(plts[[i]])
    if(i != plot.type[length(plot.type)]){
      cat(paste("plot.type =",i,"\n"))
      ans <- readline(prompt="Hit <Return> to see next plot (q to exit): ")
      if(tolower(ans)=="q") return()
    }else{
      cat("Finished plots.\n")
    }
  }
}