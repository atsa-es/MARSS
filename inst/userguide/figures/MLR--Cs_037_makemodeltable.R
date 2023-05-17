###################################################
### code chunk number 40: Cs_037_makemodeltable
###################################################
if (!exists("tabledir")) tabledir <- ""
slope.names <- paste("D", rownames(dat), sep = ".")
phi.names <- names(coef(sleep.mod4, type = "vector"))[str_detect(names(coef(sleep.mod4, type = "vector")), "B.")]

model.data <- cbind(
  c(logLik(sleep.lm2), coef(sleep.lm2)[19:36], rep(NA, nsub)),
  c(sleep.mod2$logLik, coef(sleep.mod2, type = "vector")[slope.names], rep(NA, nsub)),
  c(sleep.mod3$logLik, coef(sleep.mod3, type = "vector")[slope.names], rep(NA, nsub)),
  c(sleep.mod4$logLik, coef(sleep.mod4, type = "vector")[c(slope.names, phi.names)]),
  c(sleep.mod5$logLik, coef(sleep.mod5, type = "vector")[c(slope.names, rep("B.diag", nsub))]),
  c(logLik(sleep.mod5.gls), coef(sleep.mod5.gls)[19:36], rep(coef(sleep.mod5.gls$modelStruct[[1]], unconstrained = FALSE), nsub))
)
rownames(model.data) <- c("logLik", paste("slope", unique(sleepstudy$Subject)), paste("phi", unique(sleepstudy$Subject)))
colnames(model.data) <- c("lm", "mod2 em", "mod3 em", "mod4 em", "mod5 em", "mod5 gls")
tmpaln <- "c" # figure out the number of cols automatically
for (i in 1:ncol(model.data)) tmpaln <- paste(tmpaln, "c", sep = "")
thetable <- xtable(model.data, caption = "Parameter estimates of different versions of the model where each subject has a separate intercept (response time on normal sleep) and different slope by day (increase in response time with each day of sleep deprivation).  The model types are discussed in the text.", label = "ref:tablesleepstudy", align = tmpaln, digits = 2)
print(thetable, type = "latex", file = paste(tabledir, "tablesleepstudy.tex", sep = ""), include.rownames = TRUE, include.colnames = TRUE, caption.placement = "top", table.placement = "htp", sanitize.text.function = function(x) {
  x
}, hline.after = c(-1, 0, nrow(model.data)))


