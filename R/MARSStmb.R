MARSStmb <- function(MLEobj) {
  
out <- marssTMB::estimate_marss(MLEobj)
obj1 <- out$obj
opt1 <- out$opt
  
# Add names to the par output; par vec needs to be in this order
# In MARSS(), Dd and Cc are within the A and U matrices in marss form
# The par element of a MLEobj is in this marss form
# {MARSS} has a helper function to convert from marxss (with Dd and Cc) to marss 
parlist <- list()
# These include D, d, C, c
model.elem <- attr(MODELobj, "par.names")
for (elem in model.elem[!(model.elem %in% c("R", "Q", "V0"))]) {
  if (!MARSS:::is.fixed(MODELobj$free[[elem]])) {
    parlist[[elem]] <- matrix(opt1$par[names(opt1$par) == elem], ncol = 1)
    parnames <- paste0(elem, ".", levels(obj1$env$map[[elem]]))
    names(opt1$par)[names(opt1$par) == elem] <- parnames
    rownames(parlist[[elem]]) <- levels(obj1$env$map[[elem]])
  } else {
    parlist[[elem]] <- matrix(0, nrow = 0, ncol = 1)
  }
}
for (elem in c("R", "Q", "V0")) {
  if (!MARSS:::is.fixed(MODELobj$free[[elem]])) { # get a new par if needed
    val <- paste0("FullCovMat", elem)
    the.par <- obj1$report()[[val]]
    d <- MARSS:::sub3D(MLEobj$model$free[[elem]], t = 1)
    # A bit of a hack but I want to allow any varcov contraints (d mat)
    # Also ensures that the par names are in the right order;
    # They might not be since TMB code split out the diag separate from offdiag
    parlist[[elem]] <- solve(crossprod(d)) %*% t(d) %*% MARSS:::vec(the.par)
  } else {
    parlist[[elem]] <- matrix(0, nrow = 0, ncol = 1)
  }
}
# par is in marxss form with D, d, C, c
MLEobj[["par"]] <- parlist[model.elem]
# We need this in marss form; only.par means the marss element is already ok
MLEobj <- MARSS:::marxss_to_marss(MLEobj, only.par = TRUE)

MLEobj$iter.record <-
  list(obj.function = obj1, opt.output = opt1)
MLEobj$start <- MLEobj$start
MLEobj$convergence <- opt1$convergence

if (fun.opt == "optim" && opt1$convergence > 1) {
  if (!control$silent) cat(paste0(pkg, "() stopped with errors. No parameter estimates returned.\n"))
  MLEobj$par <- MLEobj$kf <- MLEobj$logLik <- NULL
  MLEobj$errors <- opt1$message
  return(MLEobj)
}

MLEobj$states <- obj1$env$parList()$X
rownames(MLEobj$states) <- attr(MLEobj$model, "X.names")
MLEobj$numIter <- opt1$iterations
MLEobj$logLik <- -1*opt1$objective
## Add AIC and AICc to the object
MLEobj <- MARSS::MARSSaic(MLEobj)

funinfo <- paste0("Function ", fun.opt, " used for optimization and TMB for likelihood calculation.\n")
if ((!control$silent || control$silent == 2) && opt1$convergence == 0) cat(paste0("Success! Converged in ", opt1$iterations, " iterations.\n", funinfo))
if ((!control$silent || control$silent == 2) && opt1$convergence == 1) cat(paste0("Warning! Max iterations of ", control$maxit, " reached before convergence.\n", funinfo))

return(MLEobj)
}