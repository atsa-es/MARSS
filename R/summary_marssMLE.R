###############################################################################################################################################
#  Summary method for class marssMLE. 
###############################################################################################################################################

summary.marssMLE <- function (object, digits = max(3, getOption("digits")-3), ...) 
    {
      # Call summary(marssMODEL)
      summary(object$model)

      # Call print(marssMLE)
      print(object)

      invisible(object) 
    }    
  
