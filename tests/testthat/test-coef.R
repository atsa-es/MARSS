skip_on_cran()

library(MARSS)

load("models.RData")
context("coef tests")

fits <- list()
i <- model.list$`GDP3-tinit0`
fits[[1]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, method = "BFGS")
fits[[2]] <- MARSS(i$data, model = fits[[1]], silent=TRUE)
fits[[3]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, fit=FALSE)
i <- model.list$RW_1D_b
fits[[4]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE)
fits[[5]] <- MARSS(i$data, model = fits[[4]], silent=TRUE)
fits[[6]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, fit=FALSE)
i <- model.list$`GDP1`
fits[[7]] <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, method = "BFGS")

what <- c("par", "par.se", "par.bias", "par.lowCI", "par.upCI", "start")
  
for(val in what){

  for(i in c(1:2, 5)){
    fit <- fits[[i]]
    if( val!="par"){ 
      p1 <- try(coef(MARSSparamCIs(fit), what=val))
      test_that(paste("coef class", i, val), {
        expect_true(inherits(p1, "try-error"))
      })
    }else{
      p1 <- try(coef(fit))
      test_that(paste("coef class", i, val), {
      expect_true(!inherits(p1, "try-error"))
    })
    }
  }
 
  for(i in c(3,6)){
    fit <- fits[[i]]
    if( val!="par"){ 
      p1 <- try(coef(MARSSparamCIs(fit), what=val))
    }else{
      p1 <- try(coef(fit))
    }
    test_that(paste("coef class", i, val), {
      expect_true(inherits(p1, "try-error"))
    })
  }
  
    for(i in c(4)){
      fit <- fits[[i]]
      if( !(val %in% c("par", "par.bias"))) p1 <- try(coef(MARSSparamCIs(fit), what=val))
      if (val=="par") p1 <- try(coef(fit))
      if(val=="par.bias") p1 <- try(coef(MARSSparamCIs(fit, method="parametric", nboot=2), what=val))
      test_that(paste("coef class", i, val), {
        expect_true(!inherits(p1, "try-error"))
      })
    }

  for(i in c(7)){
    fit <- fits[[i]]
    if( !(val %in% c("par", "par.bias")) ) p1 <- try(coef(MARSSparamCIs(fit, hessian.fun="optim"), what=val))
    if (val=="par") p1 <- try(coef(fit))
    if(val=="par.bias") p1 <- try(coef(MARSSparamCIs(fit, method="parametric", nboot=2), what=val))
    test_that(paste("coef class", i, val), {
      expect_true(!inherits(p1, "try-error"))
    })
  }
}
  
