for(val in what){
  cat(val)
  browser()
  if( !(val %in% c("par", "par.bias")) ) 
    p1 <- try(coef(MARSSparamCIs(fit, hessian.fun="optim"), what=val))
  if (val=="par") p1 <- try(coef(fit))
  if(val=="par.bias") p1 <- try(coef(MARSSparamCIs(fit, method="parametric", nboot=2), what=val))
}