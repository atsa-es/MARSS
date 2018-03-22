plot.marssMLE= function(x, plot.type=c("observations", "states", "model.residuals", "state.residuals"), form=c("marxss", "marss", "dfa")) {
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
    p3 = ggplot2::ggplot(df[(!is.na(df$.resids) & !is.na(df$y)),], ggplot2::aes(t, .resids)) + 
      ggplot2::geom_point(col="blue") +
      ggplot2::xlab("Time") + 
      ggplot2::ylab("Observation residuals, y - E[y]") +
      ggplot2::facet_wrap(~.rownames, scale="free_y") + 
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype=3)
    plts[["model.residuals"]] = p3
    if(identical(plot.type, "model.residuals")) return(p3)
  }  
  
  if("state.residuals" %in% plot.type) {
    # make plot of process residuals
    df = augment.marssMLE(x, "states")
    df$.rownames = paste0("State ",df$.rownames)
    p4 = ggplot2::ggplot(df[!is.na(df$.resids),], ggplot2::aes(t, .resids)) + 
      ggplot2::geom_point(col="blue") + 
      ggplot2::xlab("Time") + 
      ggplot2::ylab("State residuals, xtT - E[x]") +
      ggplot2::facet_wrap(~.rownames, scale="free_y") + 
      ggplot2::geom_hline(ggplot2::aes(yintercept=0), linetype=3)
    plts[["state.residuals"]] = p4
    if(identical(plot.type, "state.residuals")) return(p4)
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