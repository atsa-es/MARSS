############################################################################################################################
#   checkMARSSInputs()
#   This checks user inputs to MARSS().
# MARSS = function(y, inits=NULL, model=NULL, method = "kem", form = "marxss", fit=TRUE,  silent = FALSE, control = NULL )
#   checks y, duplication in inits, method, form, fit, silent
# control is checked when the is.marssMLE is called.
##########################################################################################################################
checkMARSSInputs <- function(MARSS.inputs, silent = FALSE) {
  alldefaults <- get("alldefaults", pkg_globals)

  # Check that wrapper passed in and all wrapper elements are present
  el <- c("data", "inits", "marss", "control", "method", "form", "fit", "silent", "fun.kf")
  if (!all(el %in% names(MARSS.inputs))) {
    msg <- paste(" Element", el[!(el %in% names(MARSS.inputs))], "is missing in MARSS call.\n")
    cat("\n", "Errors were caught in checkMARSSInputs \n", msg, sep = "")
    stop("Stopped in checkMARSSInputs() due to specification problem(s).\n", call. = FALSE)
  }

  # Second make sure specified method is allowed
  # The user might further restrict this in their MARSS.form() function
  # allowed.methods is speced in .onLoad
  allowed.methods <- get("allowed.methods", envir = pkg_globals)
  if (!(MARSS.inputs$method %in% allowed.methods)) {
    msg <- paste(" ", MARSS.inputs$method, "is not among the allowed methods. See ?MARSS.\n")
    cat("\n", "Errors were caught in MARSS \n", msg, sep = "")
    stop("Stopped in MARSS() due to problem(s) with required arguments.\n", call. = FALSE)
  }

  ## Alert users re 3.0+ not compatible with control$kf.x0 or $diffuse
  if (is.list(MARSS.inputs$control)) {
    if (!is.null(MARSS.inputs$control$kf.x0)) {
      stop("Stopped in checkMARSSInputs(): control$kf.x0 is no longer used in MARSS 3+.\n  Set the initial state time using the model$tinitx=1 or =0 instead.\n", call. = FALSE)
    }
    if (!is.null(MARSS.inputs$control$diffuse)) {
      stop("Stopped in checkMARSSInputs(): diffuse is no longer part of the control list in MARSS 3+.\n  Pass in diffuse in the model list instead.\n", call. = FALSE)
    }
  }

  # check that control has some default values if not passed in
  # the user might reset these in their MARSS.form function
  req.args <- c("inits", "control") # needed but not necessarily speced by user

  # the alldefaults is set in .onLoad but user might override later in MARSS.form()
  defaults <- alldefaults[[MARSS.inputs$method]] # just use the generic values in .onLoad

  ## Now set defaults if needed, first deal with case where arg not passed in all
  for (el in req.args) {
    if (is.null(MARSS.inputs[[el]])) MARSS.inputs[[el]] <- defaults[[el]]
  }

  ## Now deal with case that inits passed in is a marssMLE object
  if (inherits(MARSS.inputs[["inits"]], "marssMLE")) {
    if (is.null(MARSS.inputs[["inits"]])) {
      stop("Stopped in checkMARSSInputs() because inits must have the par element if class marssMLE.\n", call. = FALSE)
    } else {
      MARSS.inputs[["inits"]] <- coef(MARSS.inputs[["inits"]], what = "par")
      if (!is.list(MARSS.inputs[["inits"]])) stop("Stopped in checkMARSSInputs() because par element of inits$par must be a list if inits is marssMLE object.\n", call. = FALSE)
    }
  }

  ### If some elements are missing from args use the defaults
  for (el in req.args) {
    tmp <- MARSS.inputs[[el]]
    if (is.list(defaults[[el]]) && !is.list(tmp)) { # then that el must be list
      stop(paste("Stopped in checkMARSSInputs(): arg ", el, " must be passed in as a list (or left off to use defaults).\n", sep = ""), call. = FALSE)
    }
    if (!all(names(tmp) %in% names(defaults[[el]]))) {
      bad.name <- names(tmp)[!(names(tmp) %in% names(defaults[[el]]))]
      if (el == "inits") {
        extra.msg <- "For inits, you need to use a list in the same model form. Use coef(fit) to see the elements in this list.\n"
      } else {
        extra.msg <- ""
      }
      stop(paste("\nStopped in checkMARSSInputs(): elements ", bad.name, " is not allowed in arg ", el, " (misspelled?).\n", extra.msg, sep = ""), call. = FALSE)
    }
    # set defaults for any req elements that were not passed in
    passed.in <- (names(defaults[[el]]) %in% names(tmp))
    for (i in names(defaults[[el]])[!passed.in]) {
      tmp[[i]] <- defaults[[el]][[i]]
    }
    MARSS.inputs[[el]] <- tmp
  } # for over the req args; filling in missing elements

  # Start error checking
  problem <- FALSE
  msg <- c("")
  # check that data element is matrix or vector, no dataframes, and is numeric
  if (!(is.matrix(MARSS.inputs$data) || is.vector(MARSS.inputs$data) || inherits(MARSS.inputs$data, "ts"))) {
    problem <- TRUE
    msg <- c(msg, " Data must be a matrix, vector or ts/mts object (not a data frame).\n")
  }
  if (!is.numeric(MARSS.inputs$data)) {
    problem <- TRUE
    msg <- c(msg, " Data must be numeric.\n")
  }

  # check that silent is T F
  if (!(MARSS.inputs$silent %in% c(TRUE, FALSE, 2))) {
    problem <- TRUE
    msg <- c(msg, " silent must be TRUE or FALSE.\n")
  }

  # check that fit is T F
  if (!(MARSS.inputs$fit %in% c(TRUE, FALSE))) {
    problem <- TRUE
    msg <- c(msg, " fit must be TRUE or FALSE.\n")
  }

  # check that fun.kf is allowed
  if (!(MARSS.inputs$fun.kf %in% c("MARSSkfas", "MARSSkfss"))) {
    problem <- TRUE
    msg <- c(msg, " fun.kf must be MARSSkfas or MARSSkfss (in quotes).\n")
  }

  # If there are errors, then don't proceed with the rest of the checks
  if (problem) {
    cat("\n", "Errors were caught in checkMARSSInputs \n", msg, sep = "")
    stop("Stopped in checkMARSSInputs() due to specification problem(s).\n", call. = FALSE)
  }

  # check that inits list doesn't have any duplicate names
  if (any(duplicated(names(MARSS.inputs$inits)))) {
    msg <- " Some of the inits names are duplicated.\n"
    cat("\n", "Errors were caught in checkMARSSInputs \n", msg, sep = "")
    stop("Stopped in checkMARSSInputs() due to problem(s) with inits.\n", call. = FALSE)
  }

  # Warnings
  # Warn if minit set lower than min.iter.conv.test
  control <- MARSS.inputs$control
  if (!silent && is.null(control$abstol) && !is.null(control$min.iter.conv.test) && !is.null(control$minit) && is.numeric(control$min.iter.conv.test) && is.numeric(control$minit)) {
    if (control$min.iter.conv.test > control$minit) {
      warning("checkMARSSInputs: control$minit is less than control$min.iter.conv.test.\nMinimum iterations will be determined by min.iter.conv.test.", call. = FALSE)
    }
    if (control$min.iter.conv.test > control$maxit) {
      warning("checkMARSSInputs: control$maxit is less than control$min.iter.conv.test.\nNo convergence test will be computed.", call. = FALSE)
    }
  }

  MARSS.inputs
}
