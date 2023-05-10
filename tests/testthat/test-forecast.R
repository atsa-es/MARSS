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

for(i in c(3,6)){
  fit <- fits[[i]]
  p1 <- try(forecast(fit, h=10), silent=TRUE)
  test_that(paste("forecast class", i), {
    expect_true(inherits(p1, "try-error"))
  })
  p1 <- try(predict(fit), silent=TRUE)
  test_that(paste("predict class", i), {
    expect_true(inherits(p1, "try-error"))
  })
}

for(i in c(1:2,4,5,7)){
  fit <- fits[[i]]
  p1 <- try(forecast(fit, h=10), silent=TRUE)
  test_that(paste("forecast class", i), {
    expect_true(inherits(p1, "marssPredict"))
  })
  p1 <- try(predict.marssMLE(fit), silent=TRUE)
  test_that(paste("predict class", i), {
    expect_true(inherits(p1, "marssPredict"))
  })
}

# Little harder model

dat <- t(harborSealWA)
dat <- dat[2:4, ] # remove the year row

for (Q in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")) {
  for (B in c("identity", "diagonal and unequal")) {
    for (R in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")) {
      if (B == "diagonal and unequal" && (is.list(Q) || is.list(R))) next
      mod <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = dat[, 1, drop = FALSE] * 1.1)
      if (B != "identity" && R != "zero") mod$tinitx <- 1
      if (Q == "zero" && R == "zero") next
      if (Q == "zero" && B != "identity") next
      fit <- MARSS(dat, model = mod, silent = TRUE)
      p1 <- try(forecast(fit, h=10), silent=TRUE)
      test_that(paste("forecast class", Q, R, B), {
        expect_true(inherits(p1, "marssPredict"))
      })
      p1 <- try(predict(fit, newdata=list(y=t(harborSeal)[2:4,1:4])), silent=TRUE)
      test_that(paste("predict class", Q, R, B), {
        expect_true(inherits(p1, "marssPredict"))
      })
    }}}

fit <- MARSS(dat, model=list(U=array(c("u1", "u2", "u3", "u4"),dim=c(3,1,22)) ))
p1 <- predict(fit, newdata=list(t=1:4+10, y=dat[,1:4+3]))
test_that("predict1", {
  expect_true(inherits(p1, "marssPredict"))
})
p1 <- predict(fit, newdata=list(y=dat[,1:4], t=20+1:4))
test_that("predict2", {
  expect_true(inherits(p1, "marssPredict"))
})
p1 <- predict(fit, newdata=list(y=dat[,1:4], t=20+1:4, c=1:100))
test_that("predict3", {
  expect_true(inherits(p1, "marssPredict"))
})
fit <- MARSS(dat, model=list(U=array(c("u1", "u2", "u3", "u4"),dim=c(3,1,22)), c=matrix(rnorm(3*22), 3, 22) ))
p1 <- try(predict(fit, newdata=list(y=dat[,1:4], t=20+1:4)))
test_that("predict4", {
  expect_true(inherits(p1, "try-error"))
})
p1 <- try(predict(fit, newdata=list(y=dat[,1:4], t=20+1:4, c=1)))
test_that("predict5", {
  expect_true(inherits(p1, "try-error"))
})
p1 <- predict(fit, newdata=list(y=dat[,1:4], t=20+1:4, c=matrix(rnorm(1),3,4)))
test_that("predict6", {
  expect_true(inherits(p1, "marssPredict"))
})
p1 <- predict(fit, newdata=list(y=dat[,1:4], c=matrix(rnorm(1),3,4)))
test_that("predict6", {
  expect_true(inherits(p1, "marssPredict"))
})
p1 <- predict(fit, newdata=list(c=matrix(rnorm(1),3,4)))
test_that("predict7", {
  expect_true(inherits(p1, "marssPredict"))
})


