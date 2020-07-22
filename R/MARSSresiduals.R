MARSSresiduals <- function(object, ..., type=c("tT", "tt1"), normalize = FALSE, silent=FALSE) {
  type <- match.arg(type)
  if(type=="tT") val <- MARSSresiduals.tT(object, ..., normalize=normalize, silent=silent)
  if(type=="tt1") val <- MARSSresiduals.tt1(object, ..., normalize=normalize, silent=silent)
  return(val)
}
