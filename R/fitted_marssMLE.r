###############################################################################################################################################
#  fitted method for class marssMLE
##############################################################################################################################################
fitted.marssMLE <- function(object, 
    type = c("ytT", "xtT", "xtt", "xtt1", "ytt", "ytt1"),
    interval = c("confidence", "prediction", "none"),
    level = 0.95, 
    fit.type = c("prediction", "estimate"),
    ...) {
  
  type <- match.arg(type)
  interval <- match.arg(interval)
  fit.type <- match.arg(fit.type)

  if(fit.type=="prediction") ret <- MARSSpredict(object, type=type, interval=interval, level=level, ..., output="data.frame")
  if(fit.type=="estimate") ret <- MARSSest(object, type=type, interval=interval, level=level, ...)
  ret
}
