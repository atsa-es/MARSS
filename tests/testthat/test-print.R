skip_on_cran()

library(MARSS)

load("models.RData")
context("Printing tests")

fits <- list()
i <- model.list$`GDP3-tinit0`
fits[[1]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, method = "BFGS")
fits[[2]] <- MARSS(i$data, model = fits[[1]], silent=TRUE)
fits[[3]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, fit=FALSE)
i <- model.list$RW_1D_b
fits[[4]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE)
fits[[5]] <- MARSS(i$data, model = fits[[4]], silent=TRUE)
fits[[6]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, fit=FALSE)


for(i in 1:length(fits)){
  fit <- fits[[i]]
  what <- list("fit", "model", "par", "logLik", "paramvector", "states", "data", "ytT", "states.se", "model.residuals", "state.residuals", "kfs", "Ey", "states.cis", "start", names(fit$model), names(fit$model$fixed), c("par", "logLik"), c("R", "x0"))
  
for(val in what){
  p1 <- print(fit, what=val, silent=TRUE)
  test_that(paste("print fit", i, paste(val, collapse=", ")), {
    expect_true(!inherits(p1, "try-error"))
  })
}
}

for(i in c(1,4)){
  fit <- fits[[i]]
  tmp <- try(MARSSparamCIs(fit), silent=TRUE)
  if(inherits(tmp, "try-error")){
    fit <- MARSSparamCIs(fit, hessian.fun="optim")
  }else{ fit <- tmp }
  what <- list("par.se", "par.upCI", "par.lowCI")
  
  for(val in what){
    p1 <- print(fit, what=val, silent=TRUE)
    test_that(paste("print fit", i, paste(val, collapse=", ")), {
      expect_true(!inherits(p1, "try-error"))
    })
  }
}

# zz <- file("all.Rout", open = "wt")
# sink(zz)
# for(i in 1:length(fits)){
#   fit <- fits[[i]]
#   what <- c("fit", "par", "logLik", "paramvector", "par.se", "par.upCI", "par.lowCI", "states", "data", "ytT", "states.se", "model.residuals", "state.residuals", "kfs", "Ey", "states.cis", "start", names(fit$model), names(fit$model$fixed))
#   
#   for(val in what){ cat(i, val, "\n")
#     p1 <- print(fit, what=val)
#   } }
# sink()
# file.show("all.Rout")

# for (i in model.list) {
#   method <- ifelse(stringr::str_detect(i$name, "GDP"), "BFGS", "kem")
#   fit <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, method = method)
#   print(fit)
#   readline("ok?")
# }
