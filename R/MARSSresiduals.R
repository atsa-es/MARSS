MARSSresiduals <- function(object, ..., type=c("tT", "tt1", "tt"), normalize = FALSE, silent=FALSE, fun.kf = c("MARSSkfas", "MARSSkfss")) {
  type <- match.arg(type)
  fun.kf <- match.arg(fun.kf)
  if(type=="tT") val <- MARSSresiduals.tT(object, ..., normalize=normalize, silent=silent, fun.kf=fun.kf)
  if(type=="tt1") val <- MARSSresiduals.tt1(object, ..., normalize=normalize, silent=silent, fun.kf=fun.kf)
  if(type=="tt") val <- MARSSresiduals.tt(object, ..., normalize=normalize, silent=silent, fun.kf=fun.kf)
  return(val)
}
