plot.marssMLE <- 
  function(x, 
           plot.type=c("observations", "states", "model.residuals", "state.residuals", "model.residuals.qqnorm", "state.residuals.qqnorm"), 
           form=c("marxss", "marss", "dfa")) {
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
    p1 = ggplot2::ggplot(states, ggplot2::aes(t, estimate)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax=conf.high), alpha=0.3, col="grey") +
      ggplot2::geom_line() + 
      ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
      ggplot2::facet_wrap(~term, scale="free_y")
    plts[["states"]] = p1
    if(identical(plot.type, "states")) return(p1)
  }
  
  if("observations" %in% plot.type) {
    # make plot of observations
    df = augment.marssMLE(x, "observations")
    p2 = ggplot2::ggplot(df, ggplot2::aes(t, .fitted)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .fitted-1.96*.se.fit, ymax=.fitted+1.96*.se.fit), alpha=0.3, col="grey") +
      ggplot2::geom_line() + 
      ggplot2::xlab("Time") + ggplot2::ylab("Estimate") +
      ggplot2::facet_wrap(~.rownames, scale="free_y") + 
      ggplot2::geom_point(data=df[!is.na(df$y),], ggplot2::aes(t, y), col="blue")
    plts[["observations"]] = p2
    if(identical(plot.type, "observations")) return(p2)
  }
  
  if("model.residuals" %in% plot.type) {
    # make plot of observation residuals
    df = augment.marssMLE(x, "observations")
    p1 = ggplot2::ggplot(df[(!is.na(df$.resids) & !is.na(df$y)),], ggplot2::aes(t, .resids)) + 
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
    p1 = ggplot2::ggplot(df[!is.na(df$.resids),], ggplot2::aes(t, .resids)) + 
      ggplot2::geom_point(col="blue") + 
      ggplot2::stat_smooth(method="loess") +
      ggplot2::xlab("Time") + 
      ggplot2::ylab("State residuals, xtT - E[x]") +
      ggplot2::facet_wrap(~.rownames, scale="free_y") + 
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype=3)
    plts[["state.residuals"]] = p1
    if(identical(plot.type, "state.residuals")) return(p1)
  }  

  if("model.residuals.qqnorm" %in% plot.type) {
    # make plot of observation residuals
    df = augment.marssMLE(x, "observations")
    y1 <- with(df, quantile(.std.resid[!is.na(.std.resid)], c(0.25, 0.75)))
    x1 <- qnorm(c(0.25, 0.75))
    slope <- diff(y1)/diff(x1)
    int <- y1[1L] - slope * x1[1L]
    p1 = ggplot2::ggplot(df, ggplot2::aes(qqnorm(.std.resid,plot.it=FALSE)[[1]], .std.resid)) +
      ggplot2::geom_point(na.rm=TRUE) +
      ggplot2::geom_abline(slope=slope, intercept=int) +
      ggplot2::xlab("Theoretical Quantiles") + 
      ggplot2::ylab("Standardized Model Residuals") +
      ggplot2::facet_wrap(~.rownames, scale="free_y")
    plts[["model.residuals.qqnorm"]] = p1
    if(identical(plot.type, "model.residuals.qqnorm")) return(p1)
  }  
  
  if("state.residuals.qqnorm" %in% plot.type) {
    # make qqplot of state residuals
    df = augment.marssMLE(x, "states")
    df$.rownames = paste0("State ",df$.rownames)
    y1 <- with(df, quantile(.std.resid[!is.na(.std.resid)], c(0.25, 0.75)))
    x1 <- qnorm(c(0.25, 0.75))
    slope <- diff(y1)/diff(x1)
    int <- y1[1L] - slope * x1[1L]
    p1 = ggplot2::ggplot(df, ggplot2::aes(qqnorm(.std.resid, plot.it=FALSE)[[1]], .std.resid)) +
      ggplot2::geom_point(na.rm=TRUE) +
      ggplot2::geom_abline(slope=slope, intercept=int) +
      ggplot2::xlab("Theoretical Quantiles") + 
      ggplot2::ylab("Standardized State Residuals") +
      ggplot2::facet_wrap(~.rownames, scale="free_y")
    plts[["state.residuals.qqnorm"]] = p1
    if(identical(plot.type, "state.residuals.qqnorm")) return(p1)
  }    
  for(i in plot.type){
    print(plts[[i]])
    if(i != plot.type[length(plot.type)]){
      cat(paste("plot.type =",i,"\n"))
      readline(prompt="Hit <Return> to see next plot: ")
    }else{
      cat("Finished plots.\n")
    }
  }
}