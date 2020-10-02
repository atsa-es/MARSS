###################################################
### code chunk number 2: Cs00_required_libraries
###################################################
library(MARSS)
library(xtable)


###################################################
### code chunk number 3: Cs01_read_in_data
###################################################
data(lakeWAplankton)
# we want lakeWAplanktonTrans, which has been log-transformed
# and the 0s replaced with NAs
plankdat <- lakeWAplanktonTrans
years <- plankdat[, "Year"] >= 1980 & plankdat[, "Year"] < 1990
phytos <- c(
  "Cryptomonas", "Diatoms", "Greens",
  "Unicells", "Other.algae"
)
dat.spp.1980 <- plankdat[years, phytos]


###################################################
### code chunk number 4: Cs02_transpose_data
###################################################
# transpose data so time goes across columns
dat.spp.1980 <- t(dat.spp.1980)
N.ts <- nrow(dat.spp.1980)
TT <- ncol(dat.spp.1980)


###################################################
### code chunk number 5: Cs03_zscore
###################################################
Sigma <- sqrt(apply(dat.spp.1980, 1, var, na.rm = TRUE))
y.bar <- apply(dat.spp.1980, 1, mean, na.rm = TRUE)
dat.z <- (dat.spp.1980 - y.bar) * (1 / Sigma)
rownames(dat.z) <- rownames(dat.spp.1980)


###################################################
### code chunk number 6: Cs03b_zscore
###################################################
dat.z <- zscore(dat.spp.1980)


###################################################
### code chunk number 7: Cs04_plotdata
###################################################
spp <- rownames(dat.spp.1980)
par(mfcol = c(3, 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in spp) {
  plot(dat.z[i, ], xlab = "", ylab = "Abundance index", bty = "L", xaxt = "n", pch = 16, col = "blue", type = "b")
  axis(1, 12 * (0:dim(dat.spp.1980)[2]) + 1, 1980 + 0:dim(dat.spp.1980)[2])
  title(i)
}


###################################################
### code chunk number 8: Cs05_setupZ
###################################################
Z.vals <- list(
  "z11", 0, 0,
  "z21", "z22", 0,
  "z31", "z32", "z33",
  "z41", "z42", "z43",
  "z51", "z52", "z53"
)
Z <- matrix(Z.vals, nrow = N.ts, ncol = 3, byrow = TRUE)


###################################################
### code chunk number 9: Cs06_printZ
###################################################
print(Z)


###################################################
### code chunk number 10: Cs07_setupQR
###################################################
Q <- B <- diag(1, 3)


###################################################
### code chunk number 11: Cs08_setupR
###################################################
R.vals <- list(
  "r11", 0, 0, 0, 0,
  0, "r22", 0, 0, 0,
  0, 0, "r33", 0, 0,
  0, 0, 0, "r44", 0,
  0, 0, 0, 0, "r55"
)

R <- matrix(R.vals, nrow = N.ts, ncol = N.ts, byrow = TRUE)


###################################################
### code chunk number 12: Cs09_printR
###################################################
print(R)


###################################################
### code chunk number 13: Cs10_setupR_short
###################################################
R <- "diagonal and unequal"


###################################################
### code chunk number 14: Cs11_setupU
###################################################
x0 <- U <- matrix(0, nrow = 3, ncol = 1)
A <- matrix(0, nrow = 6, ncol = 1)
x0 <- U <- A <- "zero"


###################################################
### code chunk number 15: Cs12_setupx0
###################################################
V0 <- diag(5, 3)


###################################################
### code chunk number 16: Cs13_define_model_list
###################################################
dfa.model <- list(
  Z = Z, A = "zero", R = R, B = B, U = U,
  Q = Q, x0 = x0, V0 = V0
)
cntl.list <- list(maxit = 50)


###################################################
### code chunk number 17: Cs14_fit_data
###################################################
kemz.3 <- MARSS(dat.z, model = dfa.model, control = cntl.list)


###################################################
### code chunk number 20: Cs15_plotfits
###################################################
fit <- kemz.3
spp <- rownames(dat.z)
par(mfcol = c(3, 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:length(spp)) {
  plot(dat.z[i, ], xlab = "", ylab = "abundance index", bty = "L", xaxt = "n", ylim = c(-4, 3), pch = 16, col = "blue")
  axis(1, 12 * (0:dim(dat.z)[2]) + 1, 1980 + 0:dim(dat.z)[2])
  par.mat <- coef(fit, type = "matrix")
  lines(as.vector(par.mat$Z[i, , drop = FALSE] %*% fit$states + par.mat$A[i, ]), lwd = 2)
  title(spp[i])
}


###################################################
### code chunk number 21: Cs16_set_up_two_trends_echo
###################################################
model.list <- list(m = 2, R = "diagonal and unequal")
kemz.2 <- MARSS(dat.spp.1980,
  model = model.list,
  z.score = TRUE, form = "dfa", control = cntl.list
)


###################################################
### code chunk number 23: Cs17_compare_mods_2n3
###################################################
print(cbind(
  model = c("3 trends", "2 trends"),
  AICc = round(c(kemz.3$AICc, kemz.2$AICc))
),
quote = FALSE
)


###################################################
### code chunk number 25: Cs18_setupmanytrends_echo (eval = FALSE)
###################################################
## # set new control params
## cntl.list <- list(minit = 200, maxit = 5000, allow.degen = FALSE)
## # set up forms of R matrices
## levels.R <- c(
##   "diagonal and equal",
##   "diagonal and unequal",
##   "equalvarcov",
##   "unconstrained"
## )
## model.data <- data.frame(stringsAsFactors = FALSE)
## # fit lots of models & store results
## # NOTE: this will take a long time to run!
## for (R in levels.R) {
##   for (m in 1:(N.ts - 1)) {
##     dfa.model <- list(A = "zero", R = R, m = m)
##     kemz <- MARSS(dat.z,
##       model = dfa.model, control = cntl.list,
##       form = "dfa", z.score = TRUE
##     )
##     model.data <- rbind(
##       model.data,
##       data.frame(
##         R = R,
##         m = m,
##         logLik = kemz$logLik,
##         K = kemz$num.params,
##         AICc = kemz$AICc,
##         stringsAsFactors = FALSE
##       )
##     )
##     assign(paste("kemz", m, R, sep = "."), kemz)
##   } # end m loop
## } # end R loop


###################################################
### code chunk number 26: Cs19_makemodeltable
###################################################
# you must run the code to do all the models for this section
if (exists("model.data")) {
  # calculate delta-AICc
  model.data$delta.AICc <- model.data$AICc - min(model.data$AICc)
  # calculate Akaike weights
  wt <- exp(-0.5 * model.data$delta.AICc)
  model.data$Ak.wt <- wt / sum(wt)
  # sort results
  model.tbl <- model.data[order(model.data$AICc), -4]
  # drop AICc from table
  # calculate cumulative wts
  model.tbl$Ak.wt.cum <- cumsum(model.tbl$Ak.wt)
  model.tbl <- model.tbl[, -4]
}


###################################################
### code chunk number 28: Cs20_kem3_R_unconstrained
###################################################
big.maxit.cntl.list <- list(minit = 200, maxit = 5000, allow.degen = FALSE)
model.list <- list(m = 2, R = "unconstrained")
the.fit <- MARSS(dat.z, model = model.list, form = "dfa", control = big.maxit.cntl.list)


###################################################
### code chunk number 29: Cs21_varimax
###################################################
# get the inverse of the rotation matrix
Z.est <- coef(the.fit, type = "matrix")$Z
H.inv <- 1
if (ncol(Z.est) > 1) H.inv <- varimax(coef(the.fit, type = "matrix")$Z)$rotmat


###################################################
### code chunk number 30: Cs22_rotations
###################################################
# rotate factor loadings
Z.rot <- Z.est %*% H.inv
# rotate trends
trends.rot <- solve(H.inv) %*% the.fit$states


###################################################
### code chunk number 31: Cs22_rotations_cis
###################################################
# Add CIs to marssMLE object
the.fit <- MARSSparamCIs(the.fit)
# Use coef() to get the upper and lower CIs
Z.low <- coef(the.fit, type = "Z", what = "par.lowCI")
Z.up <- coef(the.fit, type = "Z", what = "par.upCI")
Z.rot.up <- Z.up %*% H.inv
Z.rot.low <- Z.low %*% H.inv
df <- data.frame(est = as.vector(Z.rot), conf.up = as.vector(Z.rot.up), conf.low = as.vector(Z.rot.low))


###################################################
### code chunk number 32: Cs23_plotfacloadings
###################################################
# plot the factor loadings
spp <- rownames(dat.z)
minZ <- 0.05
m <- dim(trends.rot)[1]
ylims <- c(-1.1 * max(abs(Z.rot)), 1.1 * max(abs(Z.rot)))
par(mfrow = c(ceiling(m / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:m) {
  plot(c(1:N.ts)[abs(Z.rot[, i]) > minZ], as.vector(Z.rot[abs(Z.rot[, i]) > minZ, i]),
    type = "h", lwd = 2, xlab = "", ylab = "", xaxt = "n", ylim = ylims, xlim = c(0, N.ts + 1)
  )
  for (j in 1:N.ts) {
    if (Z.rot[j, i] > minZ) {
      text(j, -0.05, spp[j], srt = 90, adj = 1, cex = 0.9)
    }
    if (Z.rot[j, i] < -minZ) {
      text(j, 0.05, spp[j], srt = 90, adj = 0, cex = 0.9)
    }
    abline(h = 0, lwd = 1, col = "gray")
  } # end j loop
  mtext(paste("Factor loadings on trend", i, sep = " "), side = 3, line = .5)
} # end i loop


###################################################
### code chunk number 33: Cs24_plottrends
###################################################
# get ts of trends
ts.trends <- t(trends.rot)
par(mfrow = c(ceiling(dim(ts.trends)[2] / 2), 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
# loop over each trend
for (i in 1:dim(ts.trends)[2]) {
  # set up plot area
  plot(ts.trends[, i],
    ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
    type = "n", lwd = 2, bty = "L",
    xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )
  # draw zero-line
  abline(h = 0, col = "gray")
  # plot trend line
  par(new = TRUE)
  plot(ts.trends[, i],
    ylim = c(-1.1, 1.1) * max(abs(ts.trends)),
    type = "l", lwd = 2, bty = "L",
    xlab = "", ylab = "", xaxt = "n"
  )
  # add panel labels
  mtext(paste("Trend", i, sep = " "), side = 3, line = 0.5)
  axis(1, 12 * (0:dim(dat.spp.1980)[2]) + 1, 1980 + 0:dim(dat.spp.1980)[2])
} # end i loop (trends)


###################################################
### code chunk number 34: Cs25_func_get_DFA_fits
###################################################
# If there were no missing values, this function will return the fits and CIs
getDFAfits <- function(MLEobj, alpha = 0.05, covariates = NULL) {
  fits <- list()
  Ey <- MARSShatyt(MLEobj) # for var() calcs
  ZZ <- coef(MLEobj, type = "matrix")$Z # estimated Z
  nn <- nrow(ZZ) # number of obs ts
  mm <- ncol(ZZ) # number of factors/states
  TT <- ncol(Ey$ytT) # number of time steps
  ## check for covars
  if (!is.null(covariates)) {
    DD <- coef(MLEobj, type = "matrix")$D
    cov_eff <- DD %*% covariates
  } else {
    cov_eff <- matrix(0, nn, TT)
  }
  ## model expectation
  fits$ex <- ZZ %*% MLEobj$states + cov_eff
  ## Var in model fits
  VtT <- MARSSkfss(MLEobj)$VtT
  VV <- NULL
  for (tt in 1:TT) {
    RZVZ <- coef(MLEobj, type = "matrix")$R - ZZ %*% VtT[, , tt] %*% t(ZZ)
    SS <- Ey$yxtT[, , tt] - Ey$ytT[, tt, drop = FALSE] %*% t(MLEobj$states[, tt, drop = FALSE])
    VV <- cbind(VV, diag(RZVZ + SS %*% t(ZZ) + ZZ %*% t(SS)))
  }
  SE <- sqrt(VV)
  ## upper & lower (1-alpha)% CI
  fits$up <- qnorm(1 - alpha / 2) * SE + fits$ex
  fits$lo <- qnorm(alpha / 2) * SE + fits$ex
  return(fits)
}


###################################################
### code chunk number 35: Cs25b_get_DFA_fits
###################################################
fit.b <- getDFAfits(the.fit)


###################################################
### code chunk number 36: Cs25c_plotbestfits
###################################################
# plot the fits
ylbl <- rownames(dat.z)
w_ts <- seq(dim(dat.z)[2])
par(mfcol = c(3, 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:N.ts) {
  up <- fit.b$up[i, ]
  mn <- fit.b$ex[i, ]
  lo <- fit.b$lo[i, ]
  plot(w_ts, mn,
    xlab = "", ylab = ylbl[i], xaxt = "n", type = "n", cex.lab = 1.2,
    ylim = c(min(lo), max(up))
  )
  axis(1, 12 * (0:dim(dat.z)[2]) + 1, 1980 + 0:dim(dat.z)[2])
  points(w_ts, dat.z[i, ], pch = 16, col = "blue")
  lines(w_ts, up, col = "darkgray")
  lines(w_ts, mn, col = "black", lwd = 2)
  lines(w_ts, lo, col = "darkgray")
}


###################################################
### code chunk number 37: Cs25d_plotwithaugment
###################################################
require(ggplot2)
alpha <- 0.05
theme_set(theme_bw())
d <- residuals(the.fit, type = "tT")
d$up <- qnorm(1 - alpha / 2) * d$.sigma + d$.fitted
d$lo <- qnorm(alpha / 2) * d$.sigma + d$.fitted
ggplot(data = d) +
  geom_line(aes(t, .fitted)) +
  geom_point(aes(t, value)) +
  geom_ribbon(aes(x = t, ymin = lo, ymax = up), linetype = 2, alpha = 0.2) +
  facet_grid(~.rownames) +
  xlab("Time Step") +
  ylab("Count")


###################################################
### code chunk number 38: Cs26_set-up-covar
###################################################
temp <- t(plankdat[years, "Temp", drop = FALSE])
TP <- t(plankdat[years, "TP", drop = FALSE])


###################################################
### code chunk number 39: Cs27_fit_covar_echo
###################################################
model.list <- list(m = 2, R = "unconstrained")
kemz.temp <- MARSS(dat.spp.1980,
  model = model.list, z.score = TRUE,
  form = "dfa", control = cntl.list, covariates = temp
)
kemz.TP <- MARSS(dat.spp.1980,
  model = model.list, z.score = TRUE,
  form = "dfa", control = cntl.list, covariates = TP
)
kemz.both <- MARSS(dat.spp.1980,
  model = model.list, z.score = TRUE,
  form = "dfa", control = cntl.list, covariates = rbind(temp, TP)
)


###################################################
### code chunk number 42: Cs28_covar_AICs
###################################################
print(cbind(
  model = c("no covars", "Temp", "TP", "Temp & TP"),
  AICc = round(c(
    the.fit$AICc, kemz.temp$AICc, kemz.TP$AICc,
    kemz.both$AICc
  ))
), quote = FALSE)


###################################################
### code chunk number 43: Cs29_plotbestcovarfits
###################################################
par.mat <- coef(kemz.temp, type = "matrix")
fit.b <- par.mat$Z %*% kemz.temp$states + matrix(par.mat$A, nrow = N.ts, ncol = TT)
spp <- rownames(dat.z)
par(mfcol = c(3, 2), mar = c(3, 4, 1.5, 0.5), oma = c(0.4, 1, 1, 1))
for (i in 1:length(spp)) {
  plot(dat.z[i, ], xlab = "", ylab = "abundance index", bty = "L", xaxt = "n", ylim = c(-4, 3), pch = 16, col = "blue")
  axis(1, 12 * (0:dim(dat.z)[2]) + 1, 1980 + 0:dim(dat.z)[2])
  lines(fit.b[i, ], lwd = 2)
  title(spp[i])
}


