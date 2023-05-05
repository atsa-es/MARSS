## Set up env for globals and some globals
.onLoad <- function(libname, pkgname) {
  assign("pkg_globals", new.env(), envir = parent.env(environment()))

  kem.methods <- c("kem", "EM.KFAS", "EM.KFSS")
  optim.methods <- c("BFGS", "BFGS.KFAS", "BFGS.KFSS", "BFGS.TMB")
  nlminb.methods <- c("TMB", "nlminb.TMB")
  MARSSoptim.methods <- c("BFGS", "BFGS.KFAS", "BFGS.KFSS")
  MARSStmb.methods <- c("TMB", "nlminb.TMB", "BFGS.TMB")
  # specify what function is used for what method
  allowed.methods <- c(kem.methods, optim.methods, nlminb.methods)
  # These are arguments that are required/allowed for all forms
  common.allowed.in.MARSS.call <- c("data", "inits", "control", "method", "form", "fit", "silent", "fun.kf")
  assign("kem.methods", kem.methods, pkg_globals)
  assign("optim.methods", optim.methods, pkg_globals)
  assign("nlminb.methods", nlminb.methods, pkg_globals)
  assign("MARSSoptim.methods", MARSSoptim.methods, pkg_globals)
  assign("MARSStmb.methods", MARSStmb.methods, pkg_globals)
  assign("allowed.methods", allowed.methods, pkg_globals)
  assign("common.allowed.in.MARSS.call", common.allowed.in.MARSS.call, pkg_globals)

  # Set up the generic defaults for methods and forms
  # model is required but is not here since it is form dependent so must be specified in MARSS.form function
  # MARSS.form()
  alldefaults <- list()
  ## MARSS
  alldefaults$kem <- list(
    inits = list(B = 1, U = 0, Q = 0.05, Z = 1, A = 0, R = 0.05, x0 = -99, V0 = 0.05, G = 0, H = 0, L = 0),
    control = list(
      minit = 15, maxit = 500, abstol = 0.001, trace = 0, sparse = FALSE,
      safe = FALSE, allow.degen = TRUE, min.degen.iter = 50, degen.lim = 1.0e-04,
      min.iter.conv.test = 15, conv.test.deltaT = 9, conv.test.slope.tol = 0.5,
      demean.states = FALSE
    )
  )

  alldefaults[["BFGS"]] <- 
    list(
    inits = list(B = 1, U = 0, Q = 0.05, Z = 1, A = 0, R = 0.05, x0 = -99, V0 = 0, G = 0, H = 0, L = 0),
    control = list(
      maxit = 5000, trace = 0, REPORT = NULL, reltol = NULL, fnscale = NULL,
      parscale = NULL, ndeps = NULL, alpha = NULL, beta = NULL, gamma = NULL,
      type = NULL, lmm = NULL, factr = NULL, pgtol = NULL, tmax = NULL, temp = NULL,
      lower = NULL, upper = NULL
    )
  )
  alldefaults[["BFGS.KFAS"]] <- alldefaults[["BFGS.KFSS"]] <- alldefaults[["BFGS"]]
    
  alldefaults[["TMB"]] <- 
    alldefaults[["BFGS.TMB"]] <- 
    alldefaults[["nlminb.TMB"]] <- alldefaults[["BFGS"]]
  alldefaults[["BFGS.TMB"]][["control"]][["tmb.silent"]] <- TRUE
  alldefaults[["TMB"]]$control <- # nlminb() control
    list(
    maxit = 5000, tmb.silent = TRUE, eval.max = 5000, iter.max = 5000, trace = 0, 
    abs.tol = NULL, rel.tol = NULL, x.tol = NULL, 
    xf.tol = NULL, step.min = NULL, step.max = NULL, 
    sing.tol = NULL)
  alldefaults[["nlminb.TMB"]] <- alldefaults[["TMB"]]
  
  assign("alldefaults", alldefaults, pkg_globals)
}
