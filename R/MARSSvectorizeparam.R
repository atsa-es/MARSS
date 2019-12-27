#######################################################################################################
#   MARSSvectorizeparam  function
#   Returns a vector of the ESTIMATED parameters or if vector passed in, that is put into list form for MLEobj$marss
#######################################################################################################
MARSSvectorizeparam <- function(MLEobj, parvec = NA, what = "par") {
  # This helper function  ONLY FOR marssMODEL form=marss!!
  # if parvec=NA) returns a vector version of all the estimated parameters (for use in say optim) from a mssm  model
  # if parvec is passed in) returns a marssMLE object with par fixed by parvec
  # what says what are we replacing; needs to be like par list: par, par.se, par.lowCI, etc
  MODELobj <- MLEobj[["marss"]]
  # en = c("Z","A","R","B","U","Q","x0","V0") #sets order for paramvec
  en <- attr(MODELobj, "par.names")
  free <- MODELobj[["free"]]
  fixed <- MODELobj[["fixed"]]
  param <- MLEobj[[what]] # might be NULL if not set yet
  paramvector <- NULL
  if (length(parvec) == 1 && is.na(parvec[1])) {
    for (elem in en) {
      if (dim(param[[elem]])[1] > 0) { # there are estimates
        mat.names <- colnames(free[[elem]])
        tmp <- as.vector(param[[elem]])
        mat.names <- paste(rep(elem, length(mat.names)), rep(".", length(mat.names)), mat.names, sep = "")
        names(tmp) <- mat.names
        paramvector <- c(paramvector, tmp)
      }
    }
    return(paramvector)
  } # end if parvec==NA

  else {
    parlen <- 0
    maxvec <- NULL
    par <- list()

    ## Check length(parvec) matches number of free params
    for (elem in en) {
      mx <- dim(free[[elem]])[2]
      parlen <- parlen + mx
      maxvec <- c(maxvec, mx)
    }
    if (length(parvec) != parlen) stop("Stopped in MARSSvectorizeparam(). Length of param vector does not match # of free params.\n", call. = FALSE)
    names(maxvec) <- en

    ## Fill in values, elem by elem
    for (elem in en) {
      if (dim(free[[elem]])[2] > 0) { ## if any free params
        mx <- maxvec[[elem]]
        free.names <- colnames(free[[elem]])
        ## Params for this element
        elemvec <- parvec[1:mx]
        if (is.null(names(elemvec))) {
          num.names <- free.names
        } else {
          ## remove prefix
          prefix <- paste(elem, ".", sep = "")
          num.names <- sub(prefix, "", names(elemvec))
        }
        ## check name match
        if (!all(num.names %in% free.names)) {
          stop(paste("Stopped in MARSSvectorizeparam(). parvec names don't match model$free names in parameter", elem, "\n"), call. = FALSE)
        }
        ## Remove "used" values from parvec
        parvec <- parvec[(mx + 1):length(parvec)]

        ## Match order of names in paramvec to order of names in columns of free matrix
        matchvec <- match(free.names, num.names)
        elemvec <- elemvec[matchvec] # reorder the elem vec to match colnames in free

        tmp <- matrix(elemvec, mx, 1)
        tmp <- list(tmp)
        names(tmp) <- elem
        par <- c(par, tmp)
      } # end if any free elem
      else { ## use fixed matrix
        tmp <- matrix(0, 0, 1)
        tmp <- list(tmp)
        names(tmp) <- elem
        par <- c(par, tmp)
      }
    } # end elem loop

    MLEobj[[what]] <- par
    return(MLEobj)
  } # end if parvec arg is passed in
}
