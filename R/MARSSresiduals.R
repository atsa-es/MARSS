MARSSresiduals <- function(object, ..., type=c("tT", "tt1", "tt"), normalize = FALSE, silent=FALSE, fun.kf = c("MARSSkfas", "MARSSkfss")) {
  type <- match.arg(type)
  fun.kf <- match.arg(fun.kf)
  if(type=="tt" & fun.kf=="MARSSkfas") message("fun.kf must be MARSSkfss for type='tt'. fun.kf reset to MARSSkfss.\n")
  if(type=="tT") val <- MARSSresiduals.tT(object, ..., normalize=normalize, silent=silent, fun.kf=fun.kf)
  if(type=="tt1") val <- MARSSresiduals.tt1(object, ..., normalize=normalize, silent=silent, fun.kf=fun.kf)
  if(type=="tt") val <- MARSSresiduals.tt(object, ..., normalize=normalize, silent=silent)
  return(val)
}
