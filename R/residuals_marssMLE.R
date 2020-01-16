residuals.marssMLE <- function(object, ..., conditioning=c("T", "t-1"), normalize = FALSE) {
  conditioning <- match.arg(conditioning)
  if(conditioning=="T") val <- MARSSresiduals.tT(object, ..., normalize=normalize)
  if(conditioning=="t-1") val <- MARSSresiduals.tt1(object, ..., normalize=normalize)
  return(val)
}
