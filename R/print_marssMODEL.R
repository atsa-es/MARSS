###############################################################################################################################################
#  Print MARSS model structure 
###############################################################################################################################################

print.marssMODEL <- function (x, ...) 
    {
  
      model.dims=attr(x,"model.dims")
      n = model.dims$y[1]
      m = model.dims$x[1]
      X.names=attr(x,"X.names")
      Y.names=attr(x,"Y.names")
      
      cat(paste("\nModel form is ",attr(x,"form")[1],". Model Structure is\n",
                "m: ", m," state process(es) named ",paste(X.names,collapse=" "),"\n",
                "n: ", n," observation time series named ",paste(Y.names,collapse=" "),"\n\n",sep=""))
      tmp = NULL

      ## Print model structure for each parameter
      ## This is where the form specific function call resides
      rpt = describe.marssMODEL(x)
	  
      for (el in attr(x,"par.names")) {
	# if model structure is in English, print it
  cat(el, ": ", rpt[[el]], "\n")
      }  
      invisible(rpt)
}

