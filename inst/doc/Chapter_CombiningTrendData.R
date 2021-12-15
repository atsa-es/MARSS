###################################################
### code chunk number 2: Cs000_required_libraries
###################################################
library(MARSS)


###################################################
### code chunk number 3: Cs001_readinredddata
###################################################
head(okanaganRedds)
logRedds <- log(t(okanaganRedds)[c("aerial", "ground"), ])


###################################################
### code chunk number 4: Cs002_fig1
###################################################
# Code for plotting raw Okanagan redd counts
plot(okanaganRedds[, 1], okanaganRedds[, 2],
  xlab = "", ylab = "Redd counts", main = "", col = "red", pch = 1
)
points(okanaganRedds[, 1], okanaganRedds[, 3], col = "blue", pch = 2)
legend("topleft",
  inset = 0.1, legend = c("Aerial survey", "Ground survey"),
  col = c("red", "blue"), pch = c(1, 2)
)


###################################################
### code chunk number 5: Cs003_reddmodel1
###################################################
model1 <- list()
model1$R <- "diagonal and equal"
model1$Z <- matrix(1, 2, 1)
model1$A <- "scaling"
kem1 <- MARSS(logRedds, model = model1)


###################################################
### code chunk number 6: Cs004_reddmodel2
###################################################
model2 <- model1 # model2 is based on model1
model2$R <- "diagonal and unequal"
kem2 <- MARSS(logRedds, model = model2)


###################################################
### code chunk number 7: Cs005_reddmodel3
###################################################
model3 <- list()
model3$Q <- "diagonal and equal"
model3$R <- "diagonal and equal"
model3$U <- "equal"
model3$Z <- "identity"
model3$A <- "zero"
kem3 <- MARSS(logRedds, model = model3)


###################################################
### code chunk number 8: Cs005b_aic
###################################################
c(mod1 = kem1$AICc, mod2 = kem2$AICc, mod3 = kem3$AICc)


###################################################
### code chunk number 9: Cs006_fig2
###################################################
# Code for plotting the fit from the best model
plot(okanaganRedds[, 1], logRedds[1, ], xlab = "", ylab = "Redd counts", 
     main = "", col = "red", ylim = c(0, 8)
)
points(okanaganRedds[, 1], logRedds[2, ], col = "blue", pch = 2)
lines(okanaganRedds[, 1], c(kem1$states), lty = 1, lwd = 2)
lines(okanaganRedds[, 1], c(kem1$states + 2 * kem1$states.se), lty = 1, lwd = 1, col = "grey40")
lines(okanaganRedds[, 1], c(kem1$states - 2 * kem1$states.se), lty = 1, lwd = 1, col = "grey40")


###################################################
### code chunk number 26: Cs007_birddata
###################################################
birddat <- t(kestrel[, c("British.Columbia", "Alberta", "Saskatchewan")])
head(kestrel)


###################################################
### code chunk number 27: Cs008_plot-bird-data
###################################################
# Make a plot of the three time series
plot(kestrel[, 1], kestrel[, 2], xlab = "", ylab = "Index of kestrel abundance", main = "", col = "red", ylim = c(0, 2), pch = 21)
points(kestrel[, 1], kestrel[, 3], col = "blue", pch = 22)
points(kestrel[, 1], kestrel[, 4], col = "purple", pch = 25)
legend("topright", inset = 0.1, legend = c("British Columbia", "Alberta", "Saskatchewan"), col = c("red", "blue", "purple"), pch = c(21, 22, 25))


###################################################
### code chunk number 28: Cs009_fit-bird-model-1
###################################################
model.b1=list()
model.b1$R="diagonal and equal"
model.b1$Z=matrix(1,3,1)
kem.b1 = MARSS(birddat, model=model.b1, control=list(minit=100) )


###################################################
### code chunk number 29: Cs011_fit-bird-model-2
###################################################
model.b2 <- list()
model.b2$Q <- "diagonal and equal"
model.b2$R <- "diagonal and equal"
model.b2$Z <- "identity"
model.b2$A <- "zero"
model.b2$U <- "equal"
kem.b2 <- MARSS(birddat, model = model.b2)


###################################################
### code chunk number 30: Cs013_fit-bird-model-3
###################################################
model.b3 <- model.b2 # is is based on model.b2
# all we change is the structure of Q
model.b3$Q <- "diagonal and unequal"
model.b3$U <- "unequal"
kem.b3 <- MARSS(birddat, model = model.b3)


###################################################
### code chunk number 31: Cs015_fit-bird-model-4
###################################################
model.b4 <- list()
model.b4$Q <- "diagonal and unequal"
model.b4$R <- "diagonal and equal"
model.b4$Z <- factor(c("BC", "AB-SK", "AB-SK"))
model.b4$A <- "scaling"
model.b4$U <- "unequal"
kem.b4 <- MARSS(birddat, model = model.b4)


###################################################
### code chunk number 32: Cs016_aics
###################################################
c(mod1 = kem.b1$AICc, mod2 = kem.b2$AICc, mod3 = kem.b3$AICc, mod4 = kem.b4$AICc)


###################################################
### code chunk number 33: Cs017_plot-bird-model-4-fits
###################################################
# Make a plot of the predicted trajectory, confidence intervals, and the raw data in log-space
plot(kestrel[, 1], kestrel[, 2], xlab = "", ylab = "Index of kestrel abundance", 
     main = "", col = "red", ylim = c(0, 2), pch = 21)
points(kestrel[, 1], kestrel[, 3], col = "blue", pch = 22)
points(kestrel[, 1], kestrel[, 4], col = "purple", pch = 25)
lines(kestrel[, 1], c(kem.b4$states[1, ]), lty = 3, lwd = 2, col = "red")
lines(kestrel[, 1], c(kem.b4$states[2, ]), lty = 3, lwd = 2, col = "blue")
lines(kestrel[, 1], c(kem.b4$states[2, ] + coef(kem.b4, type = "matrix")$A[3, 1]), 
      lty = 3, lwd = 2, col = "purple")
legend("topright", inset = 0.1, legend = c("British Columbia", "Alberta", "Saskatchewan"), 
       col = c("red", "blue", "purple"), pch = c(21, 22, 25))


