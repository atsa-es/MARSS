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


