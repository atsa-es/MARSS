skip_on_cran()

library(MARSS)

load("models.RData")

for(i in model.list){
  context(paste(i$name, "plotting"))
  
  method <- ifelse(stringr::str_detect(i$name, "GDP"), "BFGS", "kem")
fit <- MARSS(i$data, model=i$model, control=i$control, silent=TRUE, method=method)
if(fit$convergence==54) next
fit2 <- fit; fit2$fun.kf="MARSSkfas"
  p1 <- try(plot(fit2), silent=TRUE)
  test_that("plot fit2", {
    expect_true(!inherits(p1, "try-error"))
  })
  library(ggplot2)
  p1 <- try(autoplot(fit2), silent=TRUE)
  test_that("plot fit2", {
    expect_true(!inherits(p1, "try-error"))
  })
  
  if(!i$kfss) next
  fit1 <- fit; fit1$fun.kf="MARSSkfss"
  p1 <- try(plot(fit1), silent=TRUE)
  test_that("plot fit1", {
    expect_true(!inherits(p1, "try-error"))
  })
  p1 <- try(autoplot(fit1), silent=TRUE)
  test_that("plot fit1", {
    expect_true(!inherits(p1, "try-error"))
  })
  
  
}