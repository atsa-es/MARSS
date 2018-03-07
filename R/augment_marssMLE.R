###############################################################################################################################################
#  augment method for class marssMLE and class marssMODEL 
#  returns a data.frame that has the data (y) and inputs (c and d)
#  for a MARXSS equation
##############################################################################################################################################
augment.marssMLE <- function (x, ...) {
  model.frame(x)
}