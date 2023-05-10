skip_on_cran()

context("Test tt calculations")

library(MARSS)

load("models.RData")

for (i in model.list) {
  #context(i$name)
  if (!i$kfss) next
  method <- ifelse(stringr::str_detect(i$name, "GDP"), "BFGS", "kem")
  fit <- MARSS(i$data, model = i$model, control = i$control, silent = TRUE, method = method)
  kf1 <- MARSSkfss(fit)
  kf2 <- MARSSkfas(fit)
  list1 <- kf1[c("xtT", "VtT", "Vtt1T", "x0T", "V0T", "xtt1", "Vtt1", "xtt", "Vtt")]
  list2 <- kf2[c("xtT", "VtT", "Vtt1T", "x0T", "V0T", "xtt1", "Vtt1", "xtt", "Vtt")]
  test_that("MARSSkf", {
    expect_equal(list1, list2)
  })
  fit1 <- fit
  fit1$fun.kf <- "MARSSkfss"
  fit2 <- fit
  fit2$fun.kf <- "MARSSkfas"
  for (type in c(c("tt1", "tT", "tt"))) {
    for (stand in c("Cholesky", "marginal", "Block.Cholesky")) {
      list1 <- residuals(fit1, type = type, standardization = stand)
      list2 <- residuals(fit2, type = type, standardization = stand)
      test_that(paste("residuals", type, stand), {
        expect_equal(list1, list2)
      })
    }
  }
  for (type in c("xtT", "xtt", "xtt1", "ytT", "ytt", "ytt1")) {
    for (interval in c("none", "confidence", "prediction")) {
      if (!i$kfss && type == "xtt") next
      if (!i$kfss && type == "ytt") next
      if (type != "ytT" && interval == "prediction") next
      if (type == "ytt" && interval != "none") next
      if (type == "ytt1" && interval != "none") next
      list1 <- tsSmooth(fit1, type = type, interval = interval)
      list2 <- tsSmooth(fit2, type = type, interval = interval)
      test_that(paste("tsSmooth", type, interval), {
        expect_equal(list1, list2, tolerance = sqrt(.Machine$double.eps * 6))
      })
    }
  }
  for (type in c("ytt1", "ytT", "xtT", "ytt", "xtt1")) {
    for (interval in c("none", "confidence", "prediction")) {
      for (out in c("data.frame", "matrix")) {
        if (!i$kfss && type == "ytt") next
        if (type != "ytT" && interval == "prediction") next
        if (type == "ytt" && interval != "none") next
        if (type == "ytt1" && interval != "none") next
        list1 <- fitted(fit1, type = type, interval = interval, output = out)
        list2 <- fitted(fit2, type = type, interval = interval, output = out)
        test_that(paste("fitted", type, interval), {
          expect_equal(list1, list2, tolerance = sqrt(.Machine$double.eps * 6))
        })
      }
    }
  }
  for (type in c("ytt1", "ytT", "xtT", "ytt", "xtt1")) {
    for (interval in c("none", "confidence", "prediction")) {
      for (x0 in c("reestimate", "use.model")) {
        for (n.ahead in c(0, 10)) {
          if (!i$kfss && type == "ytt") next
          if (type != "ytT" && interval == "prediction") next
          if (type == "ytt" && interval != "none") next
          if (type == "ytt1" && interval != "none") next
          list1 <- predict(fit1, type = type, interval = interval, x0 = x0, n.ahead = n.ahead)
          list2 <- predict(fit2, type = type, interval = interval, x0 = x0, n.ahead = n.ahead)
          list1$fun.kf <- list2$fun.kf <- list1$model$fun.kf <- list2$model$fun.kf <- NULL
          test_that(paste("predict", type, interval), {
            expect_equal(list1, list2, tolerance = sqrt(.Machine$double.eps * 6))
          })
        }
      }
    }
  }
  for (type in c("ytT", "xtT", "ytt", "ytt1", "xtt", "xtt1")) {
    for (interval in c("none", "confidence", "prediction")) {
      for (x0 in c("reestimate", "use.model")) {
        for (n.ahead in c(1, 10)) {
          if (!i$kfss && type == "ytt") next
          if (type != "ytT" && interval == "prediction") next
          if (type == "ytt" && interval != "none") next
          if (type == "ytt1" && interval != "none") next
          list1 <- forecast(fit1, type = type, interval = interval, x0 = x0, h = n.ahead)
          list2 <- forecast(fit2, type = type, interval = interval, x0 = x0, h = n.ahead)
          list1$fun.kf <- list2$fun.kf <- list1$model$fun.kf <- list2$model$fun.kf <- NULL
          test_that(paste("predict", type, interval), {
            expect_equal(list1, list2, tolerance = sqrt(.Machine$double.eps * 6))
          })
        }
      }
    }
  }
}
