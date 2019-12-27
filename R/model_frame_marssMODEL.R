###############################################################################################################################################
#  model.frame method for class marssMLE and class marssMODEL
#  returns a data.frame that has the data (y) and inputs (c and d)
#  for a MARXSS equation
##############################################################################################################################################
model.frame.marssMODEL <- function(formula, ...) {
  model <- formula
  model.form <- attr(model, "form")[1]
  model.dims <- attr(model, "model.dims")
  TT <- model.dims[["x"]][2]
  f <- model[["fixed"]]
  ret <- data.frame(tid = 1:TT, t(model$data))

  # Add on covariates/inputs c and d if they exist
  if (model.form == "dfa") {
    ddims <- model.dims[["d"]]
    if (ddims[3] != 1) {
      covariates <- f[["d"]]
      dim(covariates) <- ddims[c(1, 3)]
      rownames(covariates) <- rownames(f[["d"]])
      ret <- cbind(ret, t(covariates))
    }
  }
  if (model.form == "marxss") {
    for (el in c("d", "c")) {
      if (!all(f[[toupper(el)]] == 0) || !all(f[[el]] == 0)) { # there is el
        elval <- f[[el]]
        eldims <- model.dims[[el]]
        elval <- matrix(elval, eldims[1], eldims[3])
        rownames(elval) <- rownames(f[[el]])
        ret <- cbind(ret, t(elval))
      }
    }
  }
  ret
} # end of model.frame.marssMODEL

model.frame.marssMLE <- function(formula, ...) {
  model.frame.marssMODEL(formula$model)
}
