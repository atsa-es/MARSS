skip_on_cran()

context("Comparison to StructTS")

library(MARSS)

y <- window(treering, start = 0, end = 20)

fit1 <- StructTS(y, type = "level")

# Run with fit1 estimates
vy <- var(y, na.rm = TRUE) / 100
mod.list <- list(
  x0 = matrix(y[1]), U = "zero", tinitx = 0,
  Q = matrix(fit1$coef[1]), R = matrix(fit1$coef[2]),
  V0 = matrix(1e+06 * vy)
)
fit2 <- MARSS(as.vector(y), model = mod.list, silent = TRUE)

# Now estimate the parameters
mod.list <- list(
  x0 = matrix(y[1]), U = "zero", tinitx = 0, V0 = matrix(1e+06 * vy),
  Q = matrix("s2xi"), R = matrix("s2eps")
)
fit3 <- MARSS(as.vector(y), model = mod.list, method = "BFGS", silent = TRUE)
fit4 <- MARSS(as.vector(y),
  model = mod.list,
  control = list(allow.degen = FALSE), silent = TRUE
)
fit5 <- MARSS(as.vector(y),
  model = mod.list,
  control = list(allow.degen = FALSE), fun.kf = "MARSSkfss", silent = TRUE
)

fit2$kf <- MARSSkfss(fit2)
fit3$kf <- MARSSkfss(fit3)
fit4$kf <- MARSSkfss(fit4)
fit5$kf <- MARSSkfss(fit5)
df <- data.frame(
  StructTS = fit1$fitted, fit2 = fit2$kf$xtt[1, ],
  fit.bfgs = fit3$kf$xtt[1, ], fit.em = fit4$kf$xtt[1, ],
  fit.em.kfss = fit5$kf$xtt[1, ]
)

test_that(paste("compare StructTS filter to MARSSkfss filter"), {
  expect_equal(df$level, df$fit2)
})

test_that(paste("compare BFGS and StructTS fits"), {
  expect_equal(mean(abs((df$fit.bfgs - df$level) / df$fit.bfgs)), 0.0002660084)
})

test_that(paste("compare BFGS and EM fits"), {
  expect_equal(mean(abs((df$fit.bfgs - df$fit.em) / df$fit.bfgs)), 0.01478007)
})

test_that(paste("compare EM fits"), {
  expect_equal(df$fit.em, df$fit.em.kfss)
})

# Residuals test
resids1 <- residuals(fit1)
resids2 <- residuals(fit2, type = "tt", standardization = "marginal")

test_that(paste("compare residuals with identical models"), {
  expect_equal(as.numeric(resids1), resids2$.std.resids)
})
