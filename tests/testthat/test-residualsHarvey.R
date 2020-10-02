skip_on_cran()

# test Harvey = TRUE with and without normalize

context("Residuals algorithms treering")

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

resids1 <- MARSSresiduals(fit2, type = "tT", Harvey = TRUE, fun.kf = "MARSSkfss")
resids2 <- MARSSresiduals(fit2, type = "tT", fun.kf = "MARSSkfss")
resids3 <- MARSSresiduals(fit2, type = "tT", fun.kf = "MARSSkfas")

test_that(paste("StructTS compare residuals Harvey = TRUE and FALSE"), {
  expect_equal(resids1, resids2)
})

test_that(paste("StructTS compare residuals kfss vs kfas"), {
  expect_equal(resids3, resids2)
})

resids1 <- MARSSresiduals(fit2, type = "tT", Harvey = TRUE, fun.kf = "MARSSkfss", normalize = TRUE)
resids2 <- MARSSresiduals(fit2, type = "tT", fun.kf = "MARSSkfss", normalize = TRUE)

test_that(paste("StructTS compare normalized residuals Harvey = TRUE and FALSE"), {
  expect_equal(resids1, resids2)
})

# Little harder model

context("Residuals algorithms harborseal with NAs")


dat <- t(harborSealWA)
dat <- dat[2:4, ] # remove the year row; no miss values

for (Q in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")) {
  for (B in c("identity", "diagonal and unequal")) {
    for (R in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")) {
      if (B == "diagonal and unequal" && (is.list(Q) || is.list(R))) next
      mod <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = dat[, 1, drop = FALSE] * 1.1)
      if (B != "identity" && R != "zero") mod$tinitx <- 1
      if (Q == "zero" && R == "zero") next
      if (Q == "zero" && B != "identity") next
      fit2 <- MARSS(dat[, c(1, 6:12, 14:22)], model = mod, fun.kf = "MARSSkfas", silent = TRUE)
      resids1 <- MARSSresiduals(fit2, type = "tT", Harvey = TRUE, fun.kf = "MARSSkfss")
      resids2 <- MARSSresiduals(fit2, type = "tT", fun.kf = "MARSSkfss")
      test_that(paste("harborseal test 1", Q, R, B), {
        expect_equal(resids1, resids2)
      })

      resids1 <- MARSSresiduals(fit2, type = "tT", Harvey = TRUE, fun.kf = "MARSSkfss", normalize = TRUE)
      resids2 <- MARSSresiduals(fit2, type = "tT", fun.kf = "MARSSkfss", normalize = TRUE)
      test_that(paste("harborseal test 1 normalized", Q, R, B), {
        expect_equal(resids1, resids2)
      })

      fit2 <- MARSS(dat, model = mod, fun.kf = "MARSSkfas", silent = TRUE)
      resids1 <- MARSSresiduals(fit2, type = "tT", fun.kf = "MARSSkfss")
      resids2 <- MARSSresiduals(fit2, type = "tT", fun.kf = "MARSSkfas")

      test_that(paste("harborseal test 2", Q, R, B), {
        expect_equal(resids1, resids2)
      })
    }
  }
}
