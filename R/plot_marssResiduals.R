plot.marssResiduals <-
  function(x,
           plot.type = c("all", "residuals", "qqplot", "acf"),
           conf.int = TRUE, conf.level = 0.95, decorate = TRUE,
           plot.par = list(),
           silent = FALSE) {
    tmp <- match.arg(plot.type)
    plot.type <- eval(formals(autoplot.marssMLE)$plot.type)
    plot.type <- plot.type[stringr::str_detect(plot.type, "resids")]
    ctype <- unique(x$type)
    plot.type <- plot.type[sapply(plot.type, function(x) {
      any(stringr::str_detect(x, ctype))})]
    cname <- unique(x$name)
    plot.type <- plot.type[sapply(plot.type, function(x) {
      any(stringr::str_detect(x, cname))})]
    if (tmp == "residuals") {
      plot.type <- plot.type[!sapply(plot.type, function(x) {
        any(stringr::str_detect(x, c("qqplot", "acf")))
      })]
    }
    if (tmp == "qqplot") plot.type <- plot.type[stringr::str_detect(plot.type, "qqplot")]
    if (tmp == "acf") plot.type <- plot.type[stringr::str_detect(plot.type, "acf")]

    plot.marssMLE(x, plot.type = plot.type, conf.int = conf.int, plot.par = plot.par, silent = silent)
  }
