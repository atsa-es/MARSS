###################################################
### code chunk number 39: makemodeltable
###################################################
B.names <- c("B11", "B22", "B33", "B44", "B12", "B23", "B24", "B42")
# rename kem.plank.0 and kem.plank.1 B to make keeping track of params easier
rownames(kem.plank.0$par$B) <- paste("B", rep(1:4, times = 4), rep(1:4, each = 4), sep = "")
rownames(kem.plank.1$par$B) <- paste("B", rep(1:4, times = 4), rep(1:4, each = 4), sep = "")
names.ests <- c("B11", "B22", "B33", "B44", "B12", "B23", "B24", "B42", "C11", "C21", "C32", "C42")
Ives.ests.Obs <- c(.48, .25, .74, .6, -.39, -.17, -.11, .1, .25, .25, -.14, -.045)
Ives.ests.ML <- c(0.5, .076, .77, .55, -.39, -.02, -.1, .1, .2, .32, -.13, -.048)
model.data <- cbind(
  Ives.ests.Obs,
  c(coef(kem.plank.0)$B[B.names, ], NA, NA, NA, NA),
  c(coef(kem.plank.1)$B[B.names, ], NA, NA, NA, NA),
  c(coef(kem.plank.2)$B[B.names, ], NA, NA, NA, NA),
  c(coef(kem.plank.3)$B[B.names, ], NA, NA, NA, NA),
  c(coef(kem.plank.4)$B[B.names, ], coef(kem.plank.4)$C[c("C11", "C21"), ], NA, NA),
  c(coef(kem.plank.5)$B[B.names, ], coef(kem.plank.5)$C[c("C11", "C21"), ], coef(kem.plank.5)$B[c("C32", "C42"), ])
)
rownames(model.data) <- names.ests
colnames(model.data) <- c("Ives et al.", "Model 0", "Model 1", "Model 2", "Model 3", "Model 4", "Model 5")
tmpaln <- "c" # figure out the number of cols automatically
for (i in 1:ncol(model.data)) tmpaln <- paste(tmpaln, "c", sep = "")
thetable <- xtable(model.data, caption = "The parameter estimates under the different plankton models.  Models 0 to 3 do not include covariates, so the C elements are blank.  Bij is the effect of species $i$ on species $j$. 1=large phytoplankton, 2=small phytoplankton, 3=Daphnia, 4=non-Daphnia zooplankton. The Ives et al. (2003) estimates are from their table 2 for the low planktivory lake with the observation model.", label = "ref:tableplank", align = tmpaln, digits = 2)
print(thetable, type = "latex", file = paste(tabledir, "tableplank.tex", sep = ""), include.rownames = TRUE, include.colnames = TRUE, caption.placement = "top", table.placement = "htp", sanitize.text.function = function(x) {
  x
}, hline.after = c(-1, 0, nrow(model.data)))


