plot.marssResiduals <-
  function(x,
           plot.type = c("all", "residuals", "qqplot", "acf"),
           conf.int = TRUE, conf.level = 0.95, decorate = TRUE,
           plot.par = list(),
           silent = FALSE, ...) {
    tmp <- match.arg(plot.type)
    plot.type <- eval(formals(plot.marssMLE)$plot.type)
    plot.type <- plot.type[grepl("resids", plot.type, fixed=TRUE)]
    ctype <- unique(x$type)
    plot.type <- plot.type[sapply(plot.type, function(x) {
      rev(strsplit(x, "[.]")[[1]])[1] %in% ctype})]
    cname <- unique(x$name)
    plot.type <- plot.type[sapply(plot.type, function(x) {
      any(sapply(cname, function(s) {grepl(s, x)}))})]
    if (tmp == "residuals") {
      plot.type <- plot.type[!sapply(plot.type, function(x) {
        any(grepl("qqplot", x) | grepl("acf", x))
      })]
    }
    if (tmp == "qqplot") plot.type <- plot.type[grepl("qqplot", plot.type)]
    if (tmp == "acf") plot.type <- plot.type[grepl("acf", plot.type)]

    plot.marssMLE(x, plot.type = plot.type, 
                  conf.int = conf.int, conf.level = conf.level, 
                  plot.par = plot.par, silent = silent)
  }
