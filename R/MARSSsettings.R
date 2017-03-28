kem.methods=c("kem")
optim.methods=c("BFGS","BFGS-kf")       
allowed.methods = c(kem.methods, optim.methods)
#These are arguments that are required/allowed for all forms
common.allowed.in.MARSS.call = c("data", "inits", "MCbounds", "control", "method", "form", "fit", "silent", "fun.kf")

#Set up the generic defaults for methods and forms
#model is required but is not here since it is form dependent so must be specified in MARSS.form function
#MARSS.form()
alldefaults = list()
##MARSS
alldefaults$kem = list(
      inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=0.05),
      control=list(minit=15, maxit=500, abstol=0.001, trace=0, sparse=FALSE,
                   safe=FALSE, allow.degen=TRUE, min.degen.iter=50, degen.lim=1.0e-04, MCInit=FALSE, 
                   numInits = 500, numInitSteps = 10, min.iter.conv.test=15, conv.test.deltaT=9, conv.test.slope.tol= 0.5,  demean.states=FALSE),
      MCbounds=list(B=c(0,1), U=c(-1,1), Q = c(1,1), Z=c(0,1), A=c(-1,1), R = c(1,1), x0 = c(-1,1) ) 
)

alldefaults$BFGS = alldefaults[["BFGS-kf"]]=list(
      inits=list(B=1, U=0, Q=0.05, Z=1, A=0, R=0.05, x0=-99, V0=0),
      control=list(maxit=5000, MCInit=FALSE, numInits = 500, numInitSteps = 10, trace=0,
      REPORT=NULL, reltol=NULL, fnscale=NULL, parscale=NULL, ndeps=NULL, alpha=NULL, beta=NULL, gamma=NULL, type=NULL, lmm=NULL, factr=NULL,
      pgtol=NULL, tmax=NULL, temp=NULL, lower=NULL, upper=NULL ),
      MCbounds=list(B=c(0,1), U=c(-1,1), Q = c(1,1), Z=c(0,1), A=c(-1,1), R = c(1,1), x0 = c(-1,1) ) 
)

