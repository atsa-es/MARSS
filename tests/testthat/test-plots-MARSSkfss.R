skip_on_cran()

context("MARSSkfss versus MARSSkfas plots")

# test first 13 plots.
# no fails
# testthat::snapshot_review('plots-MARSSkfss')
options(testthat.progress.max_fails = 100)

# devtools::test_active_file()

library(MARSS)
library(ggplot2)
plottypes <- eval(formals(MARSS:::autoplot.marssMLE)$plot.type)
plottypes <- plottypes[!(plottypes %in% c("all", "residuals"))]
plottypes <- c(plottypes[1:12], "residuals")

fils <- dir(file.path(here::here(), "tests", "testthat", "_snaps", "plots-MARSSkfss"), full.names = TRUE)
file.remove(fils)

# Easy model
set.seed(123)
test.num <- "easy1"
dat <- cumsum(rnorm(20))
mod.list <- list(tinitx = 1, U = "zero", R = "unconstrained", Q = "unconstrained", B = "unconstrained", x0 = matrix(dat[1] * 1.0001))
kemfit1 <- MARSS(dat, model = mod.list, silent = TRUE)
kemfit2 <- MARSS(dat, model = mod.list, fun.kf = "MARSSkfss", silent = TRUE)

for(plt in plottypes){
  conf.int <- TRUE
  if(plt %in% c("ytt", "ytt1")) conf.int <- FALSE
  # test that new image is the same
  if(exists("a1")) rm(a1)
  fil <- paste0(stringr::str_replace(plt, "[.]", "-"), "-",test.num)
  fil <- stringr::str_replace(fil, "T", "capT")
  test_that(paste(plt,"setup-simple R small tinitx=1"), {
    a1 <- autoplot(kemfit1, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
  test_that(paste(plt,"simple R small tinitx=1"), {
    a1 <- autoplot(kemfit2, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
}

test.num <- "easy2"
mod.list <- list(tinitx = 0, U = "zero", R = "unconstrained", Q = "unconstrained", B = "unconstrained")
kemfit1 <- MARSS(dat, model = mod.list, silent = TRUE)
kemfit2 <- MARSS(dat, model = mod.list, fun.kf = "MARSSkfss", silent = TRUE)

for(plt in plottypes){
  conf.int <- TRUE
  if(plt %in% c("ytt", "ytt1")) conf.int <- FALSE
  # test that new image is the same
  if(exists("a1")) rm(a1)
  fil <- paste0(stringr::str_replace(plt, "[.]", "-"), "-",test.num)
  fil <- stringr::str_replace(fil, "T", "capT")
  test_that(paste(plt,"setup-simple R small tinitx=0"), {
    a1 <- autoplot(kemfit1, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
  test_that(paste(plt,"simple R small tinitx=0"), {
    a1 <- autoplot(kemfit2, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
}

# Little harder model

dat <- t(harborSealWA)
dat <- dat[2:4, ] # remove the year row
testid <- 1
for (Q in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")) {
  for (B in c("identity", "diagonal and unequal")) {
    for (R in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")) {
      if (B == "diagonal and unequal" && (is.list(Q) || is.list(R))) next
      mod <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = dat[, 1, drop = FALSE] * 1.1)
      if (B != "identity" && R != "zero") mod$tinitx <- 1
      if (Q == "zero" && R == "zero") next
      if (Q == "zero" && B != "identity") next
      kemfit1 <- MARSS(dat, model = mod, fun.kf = "MARSSkfss", silent = TRUE)
      kemfit2 <- MARSS(dat, model = mod, fun.kf = "MARSSkfas", silent = TRUE)
      test.num <- paste0("harborseal",testid)
      testid <- testid+1
      for(plt in plottypes){
        conf.int <- TRUE
        if(plt %in% c("ytt", "ytt1")) conf.int <- FALSE
        # test that new image is the same
        if(exists("a1")) rm(a1)
        fil <- paste0(stringr::str_replace(plt, "[.]", "-"), "-",test.num)
        fil <- stringr::str_replace(fil, "T", "capT")
        test_that(paste(plt,"setup-harborseal", Q, R, B), {
          a1 <- autoplot(kemfit1, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
          vdiffr::expect_doppelganger(fil, a1)
        })
        test_that(paste(plt,"harborseal", Q, R, B), {
          a1 <- autoplot(kemfit2, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
          vdiffr::expect_doppelganger(fil, a1)
        })
      }
    }
  }
}

# test B unconstrained

Q <- "diagonal and unequal"
B <- "unconstrained"
R <- "diagonal and unequal"
mod <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = "unequal", tinitx = 1)
kemfit1 <- MARSS(dat, model = mod, fun.kf = "MARSSkfss", silent = TRUE)
kemfit2 <- MARSS(dat, model = mod, fun.kf = "MARSSkfas", silent = TRUE)
test.num <- "Bunconst1"
for(plt in plottypes){
  conf.int <- TRUE
  if(plt %in% c("ytt", "ytt1")) conf.int <- FALSE
  # test that new image is the same
  if(exists("a1")) rm(a1)
  fil <- paste0(stringr::str_replace(plt, "[.]", "-"), "-",test.num)
  fil <- stringr::str_replace(fil, "T", "capT")
  test_that(paste(plt,"setup-harborseal-Bunconst"), {
    a1 <- autoplot(kemfit1, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
  test_that(paste(plt,"harborseal-Bunconst"), {
    a1 <- autoplot(kemfit2, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
}

# test Q with some zeros

Q <- ldiag(list("q1", 0, "q2"))
B <- "identity"
R <- "diagonal and equal"
mod <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = "unequal", tinitx = 1)
kemfit1 <- MARSS(dat, model = mod, fun.kf = "MARSSkfss", silent = TRUE)
kemfit2 <- MARSS(dat, model = mod, fun.kf = "MARSSkfas", silent = TRUE)
test.num <- "Q0s"
for(plt in plottypes){
  conf.int <- TRUE
  if(plt %in% c("ytt", "ytt1")) conf.int <- FALSE
  # test that new image is the same
  if(exists("a1")) rm(a1)
  fil <- paste0(stringr::str_replace(plt, "[.]", "-"), "-",test.num)
  fil <- stringr::str_replace(fil, "T", "capT")
  test_that(paste(plt,"setup", test.num), {
    a1 <- autoplot(kemfit1, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
  test_that(paste(plt, test.num), {
    a1 <- autoplot(kemfit2, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
}

R <- ldiag(list("r1", 0, "r2"))
B <- "identity"
Q <- "diagonal and equal"
mod <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0 = "unequal", tinitx = 1)
kemfit1 <- MARSS(dat, model = mod, fun.kf = "MARSSkfss", silent = TRUE)
kemfit2 <- MARSS(dat, model = mod, fun.kf = "MARSSkfas", silent = TRUE)
test.num <- "R0s"
for(plt in plottypes){
  conf.int <- TRUE
  if(plt %in% c("ytt", "ytt1")) conf.int <- FALSE
  # test that new image is the same
  if(exists("a1")) rm(a1)
  fil <- paste0(stringr::str_replace(plt, "[.]", "-"), "-",test.num)
  fil <- stringr::str_replace(fil, "T", "capT")
  test_that(paste(plt,"setup", test.num), {
    a1 <- autoplot(kemfit1, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
  test_that(paste(plt, test.num), {
    a1 <- autoplot(kemfit2, plot.type=plt, silent=TRUE, conf.int = conf.int)[[1]]
    vdiffr::expect_doppelganger(fil, a1)
  })
}

