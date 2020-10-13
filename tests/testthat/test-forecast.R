skip_on_cran()

library(MARSS)

load("models.RData")
context("Forecast and Predict tests")

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

for(i in c(1:3,6)){
  fit <- fits[[i]]
  p1 <- try(forecast.marssMLE(fit, h=10), silent=TRUE)
  test_that(paste("forecast class", i), {
    expect_true(inherits(p1, "try-error"))
  })
  p1 <- try(predict(fit), silent=TRUE)
  test_that(paste("predict class", i), {
    expect_true(inherits(p1, "try-error"))
  })
}

for(i in c(4,5,7)){
  fit <- fits[[i]]
  p1 <- try(forecast.marssMLE(fit, h=10), silent=TRUE)
  test_that(paste("forecast class", i), {
    expect_true(inherits(p1, "marssPredict"))
  })
  p1 <- try(predict.marssMLE(fit), silent=TRUE)
  test_that(paste("predict class", i), {
    expect_true(inherits(p1, "marssPredict"))
  })
  
}

