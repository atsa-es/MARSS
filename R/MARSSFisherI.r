#######################################################################################################
#   MARSSFisherI function
#   Compute the Fisher Information Matrix via Harvey () and Cavanaugh & Shumway (1996) methods
#   Louis (1982) is nicely related to the derivatives used by the EM algorithm, but the 
#   temporally correlation must be dealt with, e.g. Duan & Fulop (2011) A stable estimator of the information matrix 
#   under EM for dependent data. Statistics and Computing 21: 83-91
#   Oakes (1999) has an approach that only involves derivatives of E(LL(theta)|data,theta') but one of the derivatives 
#   will be the derivative of the E(X|data,theta') with respect to theta'.  It's not clear how to do that derivative.
#   Ho, Shumway and Ombao (2006) The state-space approach to modeling dynamic processes in Models for Intensive Longitudinal Data
#   page 157 suggest that this derivative is hard to compute.
#######################################################################################################
MARSSFisherI = function( MLEobj, method="Harvey1989") {
  #options are Harvey1989, Optim, 
  
    if(method=="optim"){
      tmp = MARSShessian( MLEobj )
      # Q and R are TRANSFORMED here
      obsFI = tmp$Hessian
    }
  if(method=="Harvey1989"){
    obsFI = MARSSharveyobsFI( MLEobj )
  }
  
  return(obsFI)
}
