## Set up env for globals and some globals
.onLoad <- function(libname, pkgname) {
  assign("pkg_globals", new.env(), envir=parent.env(environment()))

  kem.methods=c("kem")
  optim.methods=c("BFGS","BFGS-kf")       
  allowed.methods = c(kem.methods, optim.methods)
  #These are arguments that are required/allowed for all forms
  common.allowed.in.MARSS.call = c("data", "inits", "control", "method", "form", "fit", "silent", "fun.kf")
  assign("kem.methods", kem.methods, pkg_globals)
  assign("optim.methods", optim.methods, pkg_globals)
  assign("allowed.methods", allowed.methods, pkg_globals)
  assign("common.allowed.in.MARSS.call", common.allowed.in.MARSS.call, pkg_globals)
  
  #Set up the generic defaults for methods and forms
  #model is required but is not here since it is form dependent so must be specified in MARSS.form function
  #MARSS.form()
  alldefaults = list()
  ##MARSS
  alldefaults$kem = list(
    inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=0.05, G=0, H=0, L=0),
    control=list(minit=15, maxit=500, abstol=0.001, trace=0, sparse=FALSE,
                 safe=FALSE, allow.degen=TRUE, min.degen.iter=50, degen.lim=1.0e-04, 
                 min.iter.conv.test=15, conv.test.deltaT=9, conv.test.slope.tol= 0.5,
                 demean.states=FALSE)
  )
  
  alldefaults$BFGS = alldefaults[["BFGS-kf"]]=list(
    inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=0, G=0, H=0, L=0),
    control=list(maxit=5000, trace=0, REPORT=NULL, reltol=NULL, fnscale=NULL, 
                 parscale=NULL, ndeps=NULL, alpha=NULL, beta=NULL, gamma=NULL, 
                 type=NULL, lmm=NULL, factr=NULL, pgtol=NULL, tmax=NULL, temp=NULL, 
                 lower=NULL, upper=NULL )
  )
  
  assign("alldefaults", alldefaults, pkg_globals)
}