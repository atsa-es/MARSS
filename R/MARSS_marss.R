###################################################################################
# marss form definition file
###################################################################################

###################################################################################
# MARSS.marss create marss from MARSS inputs
# marss_to_marss not needed
# print_marss
# coef_marss
# MARSSintis_marss
# predict_marss
# describe_marss
# is.marssMODEL_marss

###################################################################################
MARSS.marss <- function(MARSS.call) {
  # load need package globals
  common.allowed.in.MARSS.call <- get("common.allowed.in.MARSS.call", envir = pkg_globals)

  # Part 1 Set up defaults and check that what the user passed in is allowed

  # Check that no args were passed into MARSS that are not allowed
  marss.allowed.in.MARSS.call <- c("model")
  allowed.in.call <- c(marss.allowed.in.MARSS.call, common.allowed.in.MARSS.call)
  if (any(!(names(MARSS.call) %in% allowed.in.call))) {
    bad.names <- names(MARSS.call)[!(names(MARSS.call) %in% allowed.in.call)]
    msg <- paste("Argument ", paste(bad.names, collapse = ", "), "  not allowed MARSS call for form ", MARSS.call$form, ". See ?MARSS.marxss\n", sep = "")
    cat("\n", "Errors were caught in MARSS.marss \n", msg, sep = "")
    stop("Stopped in MARSS.marss() due to problem(s) with model specification.\n", call. = FALSE)
  }

  # 1 Check for form dependent user inputs for method and reset defaults for inits and control if desired

  # 2 Specify the text shortcuts and whether factors or matrices can be passed in
  #   The names in the allowed list do not need to be A, B, Q .... as used in form=marss object
  #   Other names can be used if you want the user to use those names; then in the MARSS.form function
  #   you convert the user passed in names into the form marss names with the A, B, Q, R, ... names
  #   checkModelList() will check what the user passes in against these allowed values, so
  #   so you need to make sure each name in model.defaults has a model.allowed value here
  model.allowed <- list(
    A = c("scaling", "unconstrained", "unequal", "equal", "zero"),
    B = c("identity", "zero", "unconstrained", "unequal", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    Q = c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    R = c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    U = c("unconstrained", "equal", "unequal", "zero"),
    x0 = c("unconstrained", "equal", "unequal", "zero", "diagonal and unequal", "diagonal and equal"),
    V0 = c("identity", "zero", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov"),
    Z = c("identity", "unconstrained", "diagonal and unequal", "diagonal and equal", "equalvarcov", "onestate"),
    G = c("identity", "zero"),
    H = c("identity", "zero"),
    L = c("identity", "zero"),
    tinitx = c(0, 1),
    diffuse = c(TRUE, FALSE),
    factors = c("Z"),
    matrices = c("A", "B", "Q", "R", "U", "x0", "Z", "V0", "G", "H", "L")
  )

  # model.defaults is form dependent so you must specify it
  model.defaults <- list(
    Z = "identity", A = "scaling", R = "diagonal and equal",
    B = "identity", U = "unconstrained", Q = "diagonal and unequal",
    x0 = "unconstrained", V0 = "zero", G = "identity", H = "identity", L = "identity",
    tinitx = 0, diffuse = FALSE
  )

  # This checks that what user passed in model list can be interpreted and converted to form marss
  # if no errors, it updates the model list by filling in missing elements with the defaults
  MARSS.call$model <- checkModelList(MARSS.call$model, model.defaults, model.allowed)

  # Part 2 Convert the model list to a marssMODEL object, form=marxss
  ## set up fixed and free elements
  fixed <- free <- list()

  model <- MARSS.call$model
  model.elem <- c("Z", "A", "R", "B", "U", "Q", "x0", "V0", "G", "H", "L")
  dat <- MARSS.call[["data"]]
  model.tsp <- attr(dat, "model.tsp")
  # Note dat is changed to matrix in MARSS()
  if (is.vector(dat)) dat <- matrix(dat, 1, length(dat))
  if(inherits(dat, "ts")){ model.tsp <- stats::tsp(dat); dat <- t(dat) }
  n <- dim(dat)[1]
  TT <- dim(dat)[2]
  if(is.null(model.tsp)) model.tsp <- c(1, TT, 1)
  if (is.null(rownames(dat))) {
    Y.names <- paste("Y", seq(1, n), sep = "") # paste(seq(1, n), sep="") #
    rownames(dat) <- Y.names
  } else {
    Y.names <- rownames(dat)
    if (any(duplicated(Y.names))) {
      for (i in Y.names[duplicated(Y.names)]) {
        loc <- (Y.names == i)
        nn <- sum(loc)
        Y.names[loc] <- paste(i, "-", 1:nn, sep = "")
        rownames(dat)[loc] <- Y.names[loc]
        MARSS.call[["data"]] <- dat
      }
    }
  }

  ## Set m based on Z specification IF Z was specified; errors will be reported later if m conflicts with other parameters
  m <- NA
  if (identical(model$Z, "unconstrained")) m <- n
  if (identical(model$Z, "equalvarcov")) m <- n
  if (identical(model$Z, "diagonal and equal")) m <- n
  if (identical(model$Z, "diagonal and unequal")) m <- n
  if (identical(model$Z, "onestate")) m <- 1
  if (identical(model$Z, "identity")) m <- n
  if (is.factor(model$Z)) m <- length(levels(model$Z))
  if (is.array(model$Z)) m <- dim(model$Z)[2]

  X.names <- NULL
  if (!(is.null(model[["X.names"]]))) X.names <- model[["X.names"]]
  if (is.null(X.names) & identical(model$Z, "identity")) {
    X.names <- paste("X.", Y.names, sep = "")
  }
  if (is.null(X.names) && is.array(model$Z)) {
    if (length(dim(model$Z)) == 3) {
      if (dim(model$Z)[3] == 1) {
        if (is.design(model$Z) & !is.null(colnames(model$Z))) {
          X.names <- colnames(model$Z)
        }
      }
    }
  }
  if (is.null(X.names) & is.factor(model$Z)) {
    X.names <- unique(as.character(model$Z))
  }
  if (is.null(X.names) & is.matrix(model$Z)) {
    if (is.design(model$Z) & !is.null(colnames(model$Z))) X.names <- colnames(model$Z)
  }
  if (is.null(X.names) & is.matrix(model$Z)) {
    if (is.identity(model$Z)) {
      X.names <- paste("X.", Y.names, sep = "")
    }
  }
  if (is.null(X.names)) X.names <- paste("X", seq(1, m), sep = "") # paste(seq(1, m),sep="")  #

  ## Set based on G and H specification
  ## error checking later will complain if conflict
  if (is.array(model$G)) g1 <- dim(model$G)[2] else g1 <- m
  if (is.array(model$H)) h1 <- dim(model$H)[2] else h1 <- n
  if (is.array(model$L)) l1 <- dim(model$L)[2] else l1 <- m

  # 3rd dim of params are set to 1 and will be reset to correct value at end
  model.dims <- list(
    data = c(n, TT), x = c(m, TT), y = c(n, TT), w = c(m, TT), v = c(n, TT),
    Z = c(n, m, 1), U = c(m, 1, 1), A = c(n, 1, 1), B = c(m, m, 1),
    Q = c(g1, g1, 1), R = c(h1, h1, 1), G = c(m, g1), H = c(n, h1), L = c(m, l1),
    x0 = c(m, 1, 1), V0 = c(l1, l1, 1)
  )

  ## Error-checking section that is specific to marss form
  # Note most error checking happens in checkMARSSInputs, checkModelList, and is.marssMLE
  # If a model elements passed in as a factors, make sure it is the correct length otherwise construction of marssMODEL object will break
  problem <- FALSE
  msg <- NULL
  ## Check model structures that are passed in as factors
  correct.factor.len <- list(Z = n, A = n, R = h1, B = m, U = m, Q = g1, x0 = m, V0 = l1)
  for (el in model.elem) {
    # if a factor then it needs to have the correct length otherwise construction of marssMODEL object will break and no NAs
    if (is.factor(model[[el]])) {
      if (length(model[[el]]) != correct.factor.len[[el]]) {
        problem <- TRUE
        msg <- c(msg, paste(" The model$", el, " is being passed as a factor and should be length ", correct.factor.len[[el]], " based on data dims. It's not. See help file.\n", sep = ""))
      }
      if (NA.in.fac <- NA %in% model[[el]]) {
        problem <- TRUE
        msg <- c(msg, paste(" NAs are not allowed in model factor for ", el, ". See help file.\n", sep = ""))
      }
    } # is factor
  } # end for (el in model.elem)

  # if el == Z then factor needs to have m levels
  if (is.factor(model$Z)) {
    if (length(levels(model$Z)) != m) {
      problem <- TRUE
      msg <- c(msg, " When Z is a factor, the number of levels must equal the number of state processes (m).\n")
    }
  }

  # Check that if A is scaling, then Z spec must lead to a design matrix
  if (identical(model$A, "scaling")) {
    if (is.array(model$Z) && length(dim(model$Z)) == 3) {
      if (dim(model$Z)[3] != 1 && !is.design(model$Z, zero.cols.ok = TRUE)) { # if it is a matrix
        problem <- TRUE
        msg <- c(msg, " If A is scaling(the default), then Z must be a time-constant design matrix:(0,1) and rowsums=1.\nYou need to specify A in your model list. You can construct a scaling A matrix and pass that in.")
      }
    }
    if (!is.array(model$Z) && !is.factor(model$Z)) { # if it is a string
      if (!(model$Z %in% c("onestate", "identity"))) {
        problem <- TRUE
        msg <- c(msg, " If A is scaling(the default), then Z must be a time-constant design matrix:(0,1) and rowsums=1.\nYou need to specify A in your model list. You can construct a scaling A matrix and pass that in.")
      }
    }
    if (is.matrix(model$Z) && !is.design(model$Z, zero.cols.ok = TRUE)) { # if it is a matrix, won't be array due to first test
      problem <- TRUE
      msg <- c(msg, " If A is scaling(the default), then Z must be a time-constant design matrix:(0,1) and rowsums=1.\nYou need to specify A in your model list. You can construct a scaling A matrix and pass that in.")
    }
  }
  if (is.array(model$x0) && length(dim(model$x0)) == 3) {
    if (dim(model$x0)[3] != 1) {
      problem <- TRUE
      msg <- c(msg, " x0 cannot be time-varying thus if x0 in model arg is 3D, the 3rd dim must equal 1.\n")
    }
  }
  if (is.array(model$V0) && length(dim(model$V0)) == 3) {
    if (dim(model$V0)[3] != 1) {
      problem <- TRUE
      msg <- c(msg, " V0 cannot be time-varying thus if V0 in model arg is 3D, the 3rd dim must equal 1.\n")
    }
  }
  if (is.array(model$L) && length(dim(model$L)) == 3) {
    if (dim(model$L)[3] != 1) {
      problem <- TRUE
      msg <- c(msg, " V0 cannot be time-varying thus if V0 in model arg is 3D, the 3rd dim must equal 1.\n")
    }
  }

  # If there are problems
  if (problem) {
    cat("\n", "Errors were caught in MARSS.marss \n", msg, sep = "")
    stop("Stopped in MARSS.marss() due to problem(s) with model specification.\n", call. = FALSE)
  }
  # end of error section


  ## Translate the text shortcuts into a marssMODEL object
  ## Translate the model structure names (shortcuts) into fixed and free
  ## fixed is a dim(1)*dim(2) X 1 vector of the fixed (intercepts) values
  ## free is a dim(1)*dim(2) X p vector of the free (betas) values for the p estimated elements

  model.elem <- c("Z", "A", "R", "B", "U", "Q", "x0", "V0", "G", "H", "L")
  if (which(model.elem == "Z") > which(model.elem == "A")) model.elem <- rev(model.elem) # Z must go first

  tmp <- list()
  for (el in model.elem) {
    tmp[[el]] <- "not assigned"
    if (el == "Z" & is.factor(model$Z)) {
      tmp[[el]] <- matrix(0, model.dims$Z[1], model.dims$Z[2])
      for (i in X.names) tmp[[el]][which(model$Z == i), which(as.vector(X.names) == i)] <- 1
    }
    if (el == "Z" & identical(model$Z, "onestate")) { # m=1
      tmp[[el]] <- matrix(1, n, 0)
    }
    if (identical(model[[el]], "identity")) {
      tmp[[el]] <- diag(1, model.dims[[el]][1])
    }
    if (identical(model[[el]], "diagonal and equal")) {
      tmp[[el]] <- array(list(0), dim = c(model.dims[[el]][1], model.dims[[el]][2]))
      diag(tmp[[el]]) <- "diag" # paste(el,"(diag)",sep="")
      if (length(tmp[[el]]) == 1) tmp[[el]][1, 1] <- el
    }
    if (identical(model[[el]], "diagonal and unequal")) {
      tmp[[el]] <- array(list(0), dim = c(model.dims[[el]][1], model.dims[[el]][2]))
      dim.mat <- model.dims[[el]][1]
      el.labs <- as.character(1:dim.mat)
      if (el %in% c("V0", "Q", "B")) el.labs <- X.names
      if (el %in% c("Z", "R")) el.labs <- Y.names
      diag(tmp[[el]]) <- paste("(", el.labs, ",", el.labs, ")", sep = "") # paste(el,"(",as.character(1:dim.mat),",",as.character(1:dim.mat),")",sep="")
      if (length(tmp[[el]]) == 1) tmp[[el]][1, 1] <- el
    }
    if (identical(model[[el]], "unconstrained") || identical(model[[el]], "unequal")) {
      tmp[[el]] <- array(NA, dim = c(model.dims[[el]][1], model.dims[[el]][2]))
      if (el %in% c("Q", "R", "V0")) { # variance-covariance matrices
        dim.mat <- model.dims[[el]][1]
        for (i in 1:dim.mat) {
          for (j in 1:dim.mat) tmp[[el]][i, j] <- tmp[[el]][j, i] <- paste("(", i, ",", j, ")", sep = "") # paste(el,"(",i,",",j,")",sep="")
        }
      } else { # not var-cov matrix
        row.name <- 1:model.dims[[el]][1]
        col.name <- 1:model.dims[[el]][2]
        if (el %in% c("U", "x0")) row.name <- X.names
        if (el == "A") row.name <- Y.names
        for (i in 1:model.dims[[el]][1]) {
          for (j in 1:model.dims[[el]][2]) {
            if (model.dims[[el]][2] > 1) {
              tmp[[el]][i, j] <- paste("(", row.name[i], ",", col.name[j], ")", sep = "")
            } # paste(el,"(",row.name[i],",",col.name[j],")",sep="")
            else {
              tmp[[el]][i, j] <- paste(row.name[i], sep = ",")
            } # paste(el,row.name[i],sep=",")
          }
        }
      }
      if (length(tmp[[el]]) == 1) tmp[[el]][1, 1] <- el
    } # unconstrained
    if (identical(model[[el]], "equalvarcov")) {
      tmp[[el]] <- array("offdiag", dim = c(model.dims[[el]][1], model.dims[[el]][2])) # array(paste(el,"(offdiag)",sep=""),dim=model.dims[[el]])
      diag(tmp[[el]]) <- "diag" # paste(el,"(diag)",sep="")
      if (length(tmp[[el]]) == 1) tmp[[el]][1, 1] <- el
    }
    if (identical(model[[el]], "equal")) {
      tmp[[el]] <- array("1", dim = c(model.dims[[el]][1], model.dims[[el]][2])) # array(el,dim=model.dims[[el]])
    }
    if (identical(model[[el]], "zero")) {
      tmp[[el]] <- array(0, dim = c(model.dims[[el]][1], model.dims[[el]][2]))
    }
    if (is.array(model[[el]])) {
      tmp[[el]] <- model[[el]]
    }
    if (el == "A" & identical(model[[el]], "scaling")) { # check above ensures that Z is design and time-invariant
      ## Construct A from fixed Z matrix
      tmp[[el]] <- matrix(list(), model.dims$A[1], model.dims$A[2])
      tmp[[el]][, 1] <- Y.names
      for (i in 1:m) {
        tmp[[el]][min(which(tmp$Z[, i] != 0)), 1] <- 0
      }
    }
    if (identical(tmp[[el]], "not assigned")) stop(paste("Stopped in MARSS.marxss(): tmp was not assigned for ", el, ".\n", sep = ""))
    free[[el]] <- convert.model.mat(tmp[[el]])$free
    fixed[[el]] <- convert.model.mat(tmp[[el]])$fixed

    # set the last dim of the model.dims since it was at a temp value to start
    model.dims[[el]][3] <- max(dim(free[[el]])[3], dim(fixed[[el]])[3])
  }


  # Set the marssMODEL form marss
  # This is the f+Dp form for the MARSS model used for user displays, printing and such
  marss_object <- list(fixed = fixed, free = free, data = dat, tinitx = model$tinitx, diffuse = model$diffuse)
  # set the attributes
  class(marss_object) <- "marssMODEL"
  attr(marss_object, "obj.elements") <- c("fixed", "free", "data", "tinitx", "diffuse")
  attr(marss_object, "form") <- "marss"
  attr(marss_object, "model.dims") <- model.dims
  attr(marss_object, "model.tsp") <- model.tsp
  # par.names are what needs to be in fixed/free pair
  attr(marss_object, "par.names") <- c("Z", "A", "R", "B", "U", "Q", "x0", "V0", "G", "H", "L")
  attr(marss_object, "X.names") <- X.names
  attr(marss_object, "Y.names") <- Y.names
  attr(marss_object, "equation") <- "x_{t}=B_{t}*x_{t-1}+U_{t}+G_{t}*w_{t}; w_{t}~MVN(0,Q_{t})\ny_{t}=Z_{t}*x_{t}+A_{t}+H_{t}*v_{t}; v_{t}~MVN(0,R_{t})"

  # Put the marss model into model and marss
  MARSS.call$model <- marss_object
  MARSS.call$marss <- marss_object

  ## Return MARSS call list with $marss and $model added
  MARSS.call
}

# the par element of a marssMLE object is already in form=marss.  No change needed if form is marss
print_marss <- function(x) {
  return(x)
}

# the par element of a marssMLE object is already in form=marss.  No change needed
coef_marss <- function(x) {
  return(x)
}

# MARSSinits needs a valid $par element with inits in marss form
# so no need to change inits
MARSSinits_marss <- function(MLEobj, inits) {
  return(inits)
}

predict_marss <- function(x, newdata, n.ahead, t.start) {
  # the predict function takes the marssMLE obj along with newdata and returns a
  # marssMODEL (form=marss) object constructed using info in newdata that is ready for use in prediction
  # x means the x (marssMLE object) for prediction

  # a marss marssMODEL object can be made into marxss, so convert to that and then use the predice_marxss function
  x.marxss <- marss_to_marxss(x, C.and.D.are.zero = TRUE)
  return(predict_marxss(x.marxss, newdata, n.ahead, t.start))
}

# describe_marss works generally with marxss form models (of which marss is one type)
describe_marss <- function(MODELobj, model.elem = NULL) {
  # This returns the structure of a model using text strings; works for form=marss or marxss
  # You can pass in model.elem to have is drop some elements of marxss (if you don't want these to print)
  fixed <- MODELobj$fixed
  free <- MODELobj$free
  model.dims <- attr(MODELobj, "model.dims")
  m <- model.dims$x[1]
  n <- model.dims$y[1]
  model.elem <- attr(MODELobj, "par.names")
  # set up the constr type list with an element for each par name and init val of "none"
  constr.type <- vector("list", length(model.elem))
  names(constr.type) <- model.elem
  constr.type[model.elem] <- "none"

  for (elem in model.elem) {
    TT.max <- max(dim(fixed[[elem]])[3], dim(free[[elem]])[3])
    if (TT.max > 1) constr.type[[elem]] <- "time-varying"
  }
  # For everything below, elem is time constant however tests are written to take a 3D matrix and it is not req that t=1

  ############ Check the other matrices
  # This part determines the form of free/fixed
  for (elem in model.elem) {
    dimm <- model.dims[[elem]][1]
    dimc <- model.dims[[elem]][2]
    while (identical(constr.type[[elem]], "none")) {
      # Z is not time varying and is a design matrix; note test above would have ensured Z is time constant
      if (elem == "Z" & is.fixed(free[[elem]]) & dim(fixed[[elem]])[3] == 1) {
        if (all(apply(fixed[[elem]], 3, is.design, dim = c(dimm, dimc)))) {
          Z.names <- c()
          X.names <- attr(MODELobj, "X.names")
          if (is.null(X.names)) X.names <- as.character(1:m)
          for (i in 1:n) Z.names <- c(Z.names, X.names[as.logical(fixed$Z[i, , 1])]) # fixed$Z will be time constant
          constr.type[[elem]] <- paste("design matrix with rows:", paste(Z.names, collapse = " "))
          break
        }
      } # break out of the while for elem
      if (is.fixed(free[[elem]])) { # it's fixed
        if (all(fixed[[elem]] == 0)) {
          constr.type[[elem]] <- "fixed and zero"
          break
        }
        if (all(fixed[[elem]] == 1)) {
          constr.type[[elem]] <- "fixed and all one"
          break
        }
        if (length(unique(fixed[[elem]])) == 1) {
          constr.type[[elem]] <- paste("fixed and all", fixed[[elem]][1])
          break
        }
        if (all(apply(fixed[[elem]], 3, is.identity, dim = c(dimm, dimc)))) {
          constr.type[[elem]] <- "identity"
        } else {
          constr.type[[elem]] <- "fixed"
        }
        break # break out of while
      }
      if (length(free[[elem]]) == 1) {
        constr.type[[elem]] <- "scalar"
        break
      }
      if (length(unique(as.vector(free[[elem]]))) == 1 & dim(free[[elem]])[2] == 1) {
        constr.type[[elem]] <- "all equal"
        break
      }
      if (all(apply(free[[elem]], 3, is.identity))) {
        constr.type[[elem]] <- "unconstrained"
        break
      }
      tmp.free <- free[[elem]]

      # The next parts are only for parameters that are not column vectors, so dim 2 !=1
      # fixed.free.to.formula requires that 3D mats have dim3=1, which will be true here since check above that time constant
      # creates a list matrix version of model
      if (model.dims[[elem]][2] != 1) {
        tmp.mat <- fixed.free.to.formula(fixed[[elem]], free[[elem]], c(dimm, dimc))
        if (elem %in% c("Q", "R", "V0")) { # variance-covariance matrices
          dimm <- sqrt(dim(fixed[[elem]])[1])
          if (
            all(fixed[[elem]] == 0) &
              all(unlist(tmp.mat[upper.tri(tmp.mat)]) == unlist(t(tmp.mat)[upper.tri(tmp.mat)])) &
              is.design(free[[elem]]) &
              length(unique(unlist(tmp.mat[upper.tri(tmp.mat)]))) == dimm * (dimm - 1) / 2 &
              length(unique(unlist(diag(tmp.mat)))) == dimm) {
            constr.type[[elem]] <- "unconstrained"
            break
          }
        }
        if (is.diagonal(tmp.mat)) {
          tmp.diag <- unique(diag(tmp.mat)) # it's a list
          tmp.diag.est <- tmp.diag[sapply(tmp.diag, is.character)]
          tmp.diag.fixed <- tmp.diag[sapply(tmp.diag, is.numeric)]
          if (length(tmp.diag.est) == dimm) {
            constr.type[[elem]] <- "diagonal and unequal"
            break
          }
          if (length(tmp.diag.est) == 1 & length(tmp.diag.fixed) == 0) {
            constr.type[[elem]] <- "diagonal and equal"
            break
          }
          if (length(tmp.diag.est) > 1 & length(tmp.diag.fixed) == 0) {
            constr.type[[elem]] <- paste("diagonal with ", length(tmp.diag.est), " groups", sep = "")
            break
          }
          if (length(tmp.diag.est) >= 1 & length(tmp.diag.fixed) > 0) {
            constr.type[[elem]] <- paste("diagonal with fixed elements and ", length(tmp.diag.est), " estimated values", sep = "")
            break
          }
        }
        if (is.equaltri(tmp.mat)) {
          if (elem %in% c("Q", "R", "V0")) { # variance-covariance matrices
            constr.type[[elem]] <- "one variance value and covariance value"
            break
          } else {
            constr.type[[elem]] <- "one diagonal value and one off-diagonal value"
            break
          }
        }
        if (is.blockdiag(tmp.mat)) {
          if (elem %in% c("Q", "R", "V0")) { # variance-covariance matrices
            constr.type[[elem]] <- "variance-covariance matrix with block diagonal structure"
            break
          } else {
            constr.type[[elem]] <- "block diagonal matrix"
            break
          }
        }
      }
      constr.type[[elem]] <- "see summary()" # not assigned to one of the above cases
    }
  } # model.elem

  for (elem in model.elem) {
    constr.type[[elem]] <- paste(constr.type[[elem]], " (", paste(model.dims[[elem]][1:2], collapse = " x "), ")", sep = "")
  }

  return(constr.type)
}

########################################################################
# is.marssMODEL_marss function
# Check that the marss object has all the parts it needs
# fixed, free, and par.names
# and that these have the proper size and form
# m is pulled from fixed$x0
########################################################################
is.marssMODEL_marss <- function(MODELobj, method = "kem") {
  msg <- NULL
  ## Set up par.names that should be in a marss model
  en <- c("Z", "A", "R", "B", "U", "Q", "x0", "V0", "G", "H", "L")

  # Check that par.names has these and only these names
  par.names <- attr(MODELobj, "par.names")
  if (!all(en %in% par.names)) {
    msg <- c(msg, "Element ", en[!(en %in% par.names)], " is missing from the par.names attribute of the model object.\n")
  }
  if (!all(par.names %in% en)) {
    msg <- c(msg, "Only ", en, "should be in the par.names attribute of the model object.\n")
  }
  model.dims <- attr(MODELobj, "model.dims")
  if (!all(en %in% names(model.dims))) {
    msg <- c(msg, "Element ", en[!(en %in% names(model.dims))], " is missing from the model.dims attribute of the model object.\n")
  }
  if (!is.null(msg)) { # rest of the tests won't work so stop now
    return(msg)
  }

  ###########################
  # Check model.dims are correct
  ###########################
  n <- dim(MODELobj$data)[1]
  TT <- dim(MODELobj$data)[2]
  m <- dim(MODELobj$fixed$x0)[1]
  g1 <- dim(MODELobj$fixed$G)[1] / m
  h1 <- dim(MODELobj$fixed$H)[1] / n
  l1 <- dim(MODELobj$fixed$L)[1] / m
  en <- c("Z", "A", "R", "B", "U", "Q", "x0", "V0", "G", "H", "L", "data", "x", "y", "w", "v")
  correct.dim1 <- c(Z = n, A = n, R = h1, B = m, U = m, Q = g1, x0 = m, V0 = l1, G = m, H = n, L = m, data = n, x = m, y = n, w = m, v = n)
  correct.dim2 <- c(Z = m, A = 1, R = h1, B = m, U = 1, Q = g1, x0 = 1, V0 = l1, G = g1, H = h1, L = l1, data = TT, x = TT, y = TT, w = TT, v = TT)
  for (elem in en) {
    ## Check for problems in the fixed/free pairs. Problems show up as TRUE
    dim.flag1 <- dim.flag2 <- FALSE

    # check dim
    dim.flag1 <- c(dim.flag1, !(model.dims[[elem]][1] == correct.dim1[[elem]]))
    dim.flag2 <- c(dim.flag2, !(model.dims[[elem]][2] == correct.dim2[[elem]]))
  }
  if (any(c(dim.flag1, dim.flag2))) { # There's a problem
    if (any(dim.flag1)) {
      msg <- c(msg, paste("Dim 1 of ", en[dim.flag1], "is incorrect. Dims should be ", correct.dim1[dim.flag1], ", for a marss model.\n"))
    }
    if (any(dim.flag2)) {
      msg <- c(msg, paste("Dim 2 of ", en[dim.flag2], "is incorrect. Dims should be ", correct.dim2[dim.flag2], ", for a marss model.\n"))
    }
    msg <- c("\nErrors were caught in is.marssMODEL_marss()\n", msg)
    return(msg)
  }

  fixed <- MODELobj$fixed
  free <- MODELobj$free

  ###########################
  # Check that x0, V0 and L are not time-varying
  ###########################
  en <- c("x0", "V0", "L")
  time.var <- NULL
  for (elem in en) {
    time.var.flag <- FALSE
    time.var.flag <- dim(fixed[[elem]])[3] != 1 || dim(free[[elem]])[3] != 1
    time.var <- c(time.var, time.var.flag)
  }
  if (any(time.var)) { # There's a problem
    msg <- c(msg, paste(en[time.var], "cannot be time-varying.  3rd dim of fixed and free must equal 1.\n"))
    msg <- c("\nErrors were caught in is.marssMODEL_marss()\n", msg)
    return(msg)
  }

  ###########################
  # Check that none of the var-cov matrices have negative values on the diagonal
  # and that there are no f+Dq elements only f+0q or 0+Dq
  # and D must be a design matrix, so no beta_1*q1 + beta_2*q2 elements
  ###########################
  en <- c("R", "Q", "V0")
  neg <- bad.var <- not.design <- NULL
  for (elem in en) {
    neg.flag <- bad.var.flag <- not.design.flag <- FALSE
    for (i in 1:max(dim(free[[elem]])[3], dim(fixed[[elem]])[3])) {
      if (dim(fixed[[elem]])[3] == 1) {
        i1 <- 1
      } else {
        i1 <- i
      }
      if (dim(free[[elem]])[3] == 1) {
        i2 <- 1
      } else {
        i2 <- i
      }
      if (is.fixed(free[[elem]][, , min(i, dim(free[[elem]])[3]), drop = FALSE])) { # this works on 3D mats
        zero.free.rows <- matrix(TRUE, correct.dim1[[elem]] * correct.dim2[[elem]], 1)
      } else {
        zero.free.rows <- apply(free[[elem]][, , i2, drop = FALSE] == 0, 1, all) # works on 3D mat
        # the requirement is that each estimated element (in p) appears only in one place in the varcov mat, but fixed rows (0 rows) are ok
        not.design.flag <- !is.design(free[[elem]][, , i2, drop = FALSE], strict = FALSE, zero.rows.ok = TRUE, zero.cols.ok = TRUE) # works on 3D if dim3=1
      }
      zero.fixed.rows <- apply(fixed[[elem]][, , i1, drop = FALSE] == 0, 1, all) # works on 3D
      fixed.mat <- unvec(fixed[[elem]][, , i1], dim = c(correct.dim1[[elem]], correct.dim2[[elem]]))
      if (any(!zero.fixed.rows & !zero.free.rows)) bad.var.flag <- TRUE # no f+Dq rows
      if (any(takediag(fixed.mat) < 0, na.rm = TRUE)) neg.flag <- TRUE # no negative diagonals
    } # end the for loop over time
    not.design <- c(not.design, not.design.flag)
    neg <- c(neg, neg.flag)
    bad.var <- c(bad.var, bad.var.flag)
  } # enf the for loop over elem
  if (any(neg)) {
    msg <- c(msg, paste("Negative values are on the diagonal of ", en[neg], ". Neg values are illegal on the diag of a var-cov matrix.\n", sep = ""))
  }
  if (any(bad.var)) {
    msg <- c(msg, paste("Fixed and estimated values are combined in some elements of ", en[bad.var], ". This is not allowed.\n", sep = ""))
  }
  if (any(not.design)) {
    msg <- c(msg, paste("The D matrices of ", en[not.design], " must be design matrices.\n", sep = ""))
  }

  ###########################
  # Check that V0, Q and R matrices are symmetric and positive-definite
  ###########################
  en <- c("R", "Q", "V0")
  pos <- symm <- NULL
  for (elem in en) {
    varcov.flag <- TRUE
    varcov.msg <- ""
    var.dim <- c(correct.dim1[[elem]], correct.dim2[[elem]])
    for (i in 1:model.dims[[elem]][3]) {
      if (dim(fixed[[elem]])[3] == 1) {
        i1 <- 1
      } else {
        i1 <- i
      }
      if (dim(free[[elem]])[3] == 1) {
        i2 <- 1
      } else {
        i2 <- i
      }
      # works on 3D if dim3=1
      par.as.list <- fixed.free.to.formula(fixed[[elem]][, , i1, drop = FALSE], free[[elem]][, , i2, drop = FALSE], var.dim) # coverts the fixed,free pair to a list matrix
      tmp <- is.validvarcov(par.as.list, method = method)
      varcov.flag <- varcov.flag & tmp$ok
      if (!tmp$ok) varcov.msg <- c(varcov.msg, paste(" ", tmp$error, "at t=", i, "\n", sep = ""))

      if (!varcov.flag) msg <- c(msg, paste("The variance-covariance matrix ", elem, " is not properly constrained.\n", sep = ""), varcov.msg)
    } # end for loop over time
  } # end for loop over elements

  ###########################
  # Check that crossprod(G), crossprod(H), crossprod(L) are invertible
  ###########################
  en <- c("G", "H", "L")
  pos <- symm <- NULL
  for (elem in en) {
    varcov.flag <- TRUE
    varcov.msg <- ""
    var.dim <- c(correct.dim1[[elem]], correct.dim2[[elem]])
    for (i in 1:model.dims[[elem]][3]) {
      if (dim(fixed[[elem]])[3] == 1) {
        i1 <- 1
      } else {
        i1 <- i
      }
      if (dim(free[[elem]])[3] == 1) {
        i2 <- 1
      } else {
        i2 <- i
      }
      # works on 3D if dim3=1
      # since G, H, and L are numeric, par.as.list will be a numeric matrix not list
      par.as.list <- fixed.free.to.formula(fixed[[elem]][, , i1, drop = FALSE], free[[elem]][, , i2, drop = FALSE], var.dim) # coverts the fixed,free pair to a list matrix
      # this requirement is mention in 4.4 in EM Derivation
      # simple test for invertibility via condition number
      condition.limit <- 1E10
      tmp <- kappa(crossprod(as.numeric(par.as.list))) < condition.limit # TRUE is good
      varcov.flag <- varcov.flag & tmp
      if (!tmp) varcov.msg <- c(varcov.msg, paste(" ", tmp$error, "at t=", i, "\n", sep = ""))

      if (!varcov.flag) msg <- c(msg, paste("The matrix t(", elem, ")%*%", elem, " must be invertible.\n", sep = ""), varcov.msg)
    } # end for loop over time
  } # end for loop over elements


  if (length(msg) == 0) {
    return(NULL)
  } else {
    msg <- c("\nErrors were caught in is.marssMODEL_marss()\n", msg)
    return(msg)
  }
}

