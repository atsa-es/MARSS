skip_on_cran()

library(MARSS)

# Easy model
context("Comparison to kfas small R")
set.seed(123)
dat <- cumsum(rnorm(20))
kemfit1 <- MARSS(dat, model=list(tinitx=1, U="zero", R="unconstrained", Q="unconstrained", B="unconstrained", x0=matrix(dat[1]*1.0001)), silent=TRUE)
kemfit2 <- MARSS(dat, model=list(tinitx=1, U="zero", R="unconstrained", Q="unconstrained", B="unconstrained", x0=matrix(dat[1]*1.0001)), fun.kf="MARSSkfss", silent=TRUE)

test_that("compare logLik simple R small tinitx=1", {
  expect_equal(kemfit1$logLik, kemfit2$logLik)
})

kemfit1 <- MARSS(dat, model=list(tinitx=0, U="zero", R="unconstrained", Q="unconstrained", B="unconstrained"), silent=TRUE)
kemfit2 <- MARSS(dat, model=list(tinitx=0, U="zero", R="unconstrained", Q="unconstrained", B="unconstrained"), fun.kf="MARSSkfss", silent=TRUE)

test_that("compare logLik simple R small tinitx=0", {
  expect_equal(kemfit1$logLik, kemfit2$logLik)
})

kf1 <- MARSSkfss(kemfit1)
kf2 <- MARSSkfas(kemfit1)
kf1_list <- kf1[c("xtT", "VtT", "Vtt1T", "x0T", "V0T", "xtt1", "Vtt1")]
kf2_list <- kf2[c("xtT", "VtT", "Vtt1T", "x0T", "V0T", "xtt1", "Vtt1")]
  test_that("compare kf list simple R small", {
    expect_equal(kf1_list, kf2_list)
  })


# Little harder model
context("Comparison to kfas harborSeal")

dat <- t(harborSealWA)
dat <- dat[2:4, ] # remove the year row

for( Q in list("unconstrained", "diagonal and equal", "equalvarcov", "zero"))
  for( B in c("identity", "diagonal and unequal"))
    for( R in list("unconstrained", "diagonal and equal", "equalvarcov", "zero")){
      if(B=="diagonal and unequal" && (is.list(Q) || is.list(R))) next
      mod <- list(Q = Q, Z = "identity", R = R, B = B, U = "zero", x0=dat[,1,drop=FALSE]*1.1)
      if ( B != "identity" && R != "zero") mod$tinitx = 1
      if ( Q == "zero" && R == "zero") next
      if ( Q == "zero" && B != "identity") next
      kemfit1 <- MARSS(dat, model = mod, fun.kf = "MARSSkfss", silent=TRUE)
      kemfit2 <- MARSS(dat, model = mod, fun.kf = "MARSSkfas", silent=TRUE)

      test_that(paste("compare logLik", Q, R, B), {
  expect_equal(kemfit1$logLik, kemfit2$logLik)
})
      kf1 <- MARSSkfss(kemfit1)
      kf2 <- MARSSkfas(kemfit1)
      kf1_list <- kf1[c("xtT", "VtT", "Vtt1T", "x0T", "V0T", "xtt1", "Vtt1", "logLik")]
      kf2_list <- kf2[c("xtT", "VtT", "Vtt1T", "x0T", "V0T", "xtt1", "Vtt1", "logLik")]
      test_that(paste("compare kfss to kfas", Q, R, B), {
         expect_equal(kf1_list, kf2_list)
       })
    }


# Wonky model
context("Comparison to kfas Structural")

# 1) Define some data

df_marss <- matrix(NA, 2, 10)
df_marss[1,] <- c(NA, NA, NA, -0.002666915, NA, NA, -0.002064963, NA, NA, 0.01564208)
df_marss[2,] <- c(NA, 0.0005053405, 0.001147921, -0.002476667, 0.003195476, 0.003941519, -0.001529331, 0.004960794, 0.005527753, 0.004705563)

# 2) Define State Space matrices

# Matrix Z

Z <- matrix(list("0.33*z1","z2",
                 "0.67*z1",0,
                 "z1",0,
                 "0.67*z1",0,
                 "0.33*z1",0,
                 1/3,0,
                 2/3,0,
                 1,0,
                 2/3,0,
                 1/3,0,
                 0,1,
                 0,0),2,12)

m <- nrow(Z)
p <- ncol(Z)

# Matrix R

R <- matrix(list(0),m,m)

# Matrix B

B <- matrix (list("b1",1,0,0,0,0,0,0,0,0,0,0,
                  "b2",0,1,0,0,0,0,0,0,0,0,0,
                  0,0,0,1,0,0,0,0,0,0,0,0,
                  0,0,0,0,1,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"b6",1,0,0,0,0,0,
                  0,0,0,0,0,"b7",0,1,0,0,0,0,
                  0,0,0,0,0,0,0,0,1,0,0,0,
                  0,0,0,0,0,0,0,0,0,1,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"b11",1,
                  0,0,0,0,0,0,0,0,0,0,"b12",0),12,12) 


# Matrix Q

Q <- matrix (list("q1",0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,"q6",0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,0,0,
                  0,0,0,0,0,0,0,0,0,0,"q11",0,
                  0,0,0,0,0,0,0,0,0,0,0,0),12,12)

# Rest of matrices

x0 <-  matrix(0,p,1)
A  <- matrix(0,m,1)
U <- matrix(0,p,1)
V0 <- 5*diag(1,p)
U <-  matrix(0,p,1)

# 3) Estimation

# Define model

model.gen <- list(Z=Z,A=A,R=R,B=B,U=U,Q=Q,x0=x0,V0=V0,tinitx=0)

fit <- MARSS(df_marss, model=model.gen, method="BFGS", fun.kf="MARSSkfas")
kf1 <- MARSSkfss(kemfit1, smoother=FALSE)
kf2 <- MARSSkfas(kemfit1)
kf1_list <- kf1[c("xtt1", "Vtt1", "logLik")]
kf2_list <- kf2[c("xtt1", "Vtt1", "logLik")]
test_that(paste("compare kfss to kfas structural"), {
  expect_equal(kf1_list, kf2_list) })
  