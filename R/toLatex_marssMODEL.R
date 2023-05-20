############################################
# Add to NAMESPACE
# S3method(toLatex, marssMODEL)
# S3method(toLatex, marssMLE)
#############################################
toLatex.marssMLE <- function(object, ..., file = NULL, digits = 2, greek = TRUE, orientation = "landscape", math.sty = "amsmath", output = c("pdf", "tex", "rawtex"), replace = TRUE, simplify = TRUE) {
  toLatex.marssMODEL(object$model, ..., file = file, digits = digits, greek = greek, orientation = orientation, math.sty = math.sty, output = output, replace = replace, simplify = simplify)
}
toLatex.marssMODEL <- function(object, ..., file = NULL, digits = 2, greek = TRUE, orientation = "landscape", math.sty = "amsmath", output = c("pdf", "tex", "rawtex"), replace = TRUE, simplify = TRUE) {
  output <- match.arg(output)
  # by default uses geometry package and amsmath
  requireNamespace("Hmisc")
  x <- object
  # set the environment of the subfunction to this function since I want them to have access to the passed in variables
  # environment(build.eqn.tex)=environment()
  # environment(get.mat.tex)=environment()
  # environment(parameters.time.varying)=environment()

  # set up the headers and footers
  eqn.head <- "\\begin{gather*}"
  eqn.foot <- "\\end{gather*}"
  mat.head <- "\\begin{bmatrix}"
  mat.foot <- "\\end{bmatrix}"

  marssMODEL <- x # used below
  # get the model from x (x is a marssMODEL object)
  fixed <- x$fixed
  free <- x$free
  # btw, this is not robust to name changes
  m <- dim(fixed$x0)[1]
  n <- dim(x$data)[1]
  TT <- dim(x$data)[2]
  model.dims <- attr(x, "model.dims")
  model.el <- attr(x, "par.names")

  # read in the equation from the attributes.  This will be used to construct the tex
  # needs to look like x=A*x+B*C*D+E+..+w, w~MVN(E+F+..,G+H+...)\n next eqn
  eqns <- attr(x, "equation")
  eqns <- strsplit(eqns, "\n")[[1]]

  # Build a vector of all the special names in all the equations; lhs of = and ~
  # Do some testing too and determine if time will be expanded
  eqn.special <- c()
  # the default behavior is no time expansion
  eqn.ts <- 1
  time.expand <- FALSE
  for (i in 1:length(eqns)) {
    eqn <- eqns[i]
    if (!str_detect(eqn, "=")) stop("toLatex.marssMODEL: the atttibute, equation, must have the form x=... or y=....")
    eqn.name <- str_trim(strsplit(eqn, "=")[[1]][1])
    # the part to the left of , is the deterministic part + random bit in last term
    # the eqn for the errors is after the ,
    eqn.stoc <- str_trim(strsplit(eqn, ";")[[1]][2])
    # the name is before the ~
    err.name <- str_trim(strsplit(eqn.stoc, "~")[[1]][1])
    # the eqn.name (like x or y) and the err.name (like w or v) are treated differently
    eqn.special <- c(eqn.special, str_trim(strsplit(eqn.name, "_")[[1]][1]), str_trim(strsplit(err.name, "_")[[1]][1]))

    if (!str_detect(eqn, "=")) stop("toLatex.marssMODEL: the atttibute, equation, must have the form x=... or y=....")
    eqn.name <- str_trim(strsplit(eqn, "=")[[1]][1])
    # the part to the left of , is the deterministic part + random bit in last term
    eqn.det <- str_trim(strsplit(eqn, ";")[[1]][1])
    # the eqn for the errors is after the ,
    eqn.stoc <- str_trim(strsplit(eqn, ";")[[1]][2])
    # the name is before the ~
    err.name <- str_trim(strsplit(eqn.stoc, "~")[[1]][1])
    # check that the name of the errors matches before and after the ,
    tmp <- strsplit(eqn.det, "[+]")[[1]] # tmp[length(tmp)] should be error part
    tmp <- strsplit(tmp[length(tmp)], "[*]")[[1]] # tmp[length(tmp)] should be err.name
    if (err.name != str_trim(tmp[length(tmp)])) stop("toLatex.marssMODEL: in the attribute, equation, the name of the errors after the , does not match the name before the ,")
    # these 2 lines get the bit between ( and ,
    eqn.stoc.mean <- strsplit(eqn.stoc, ",")[[1]][1]
    eqn.stoc.mean <- str_trim(strsplit(eqn.stoc.mean, "[(]")[[1]][2])
    # these 2 lines get the bit between , and )
    eqn.stoc.var <- strsplit(eqn.stoc, "[)]")[[1]][1]
    eqn.stoc.var <- str_trim(strsplit(eqn.stoc.var, ",")[[1]][2])
    # If any of the equation is time-varying then will have 1 equation for each time
    # Determine if that is necessary
    if (!simplify && (parameters.time.varying(eqn.det, eqn.special, marssMODEL) |
      parameters.time.varying(eqn.stoc.mean, eqn.special, marssMODEL) |
      parameters.time.varying(eqn.stoc.var, eqn.special, marssMODEL))) {
      eqn.ts <- TT
      time.expand <- TRUE
    }
  }

  all.t.eqn.tex <- c() # this will be all the equations
  for (t in 1:eqn.ts) { # one eqn for each time step; default is to simplify and not time expand
    eqn.tex <- c() # this is just for time=t
    # eqns is from the equation attribute of marssMODEL
    for (i in 1:length(eqns)) {
      eqn <- eqns[i]
      ##############################################################
      # First read in and parse the equation attribute of marssMODEL
      ##############################################################
      eqn.name <- str_trim(strsplit(eqn, "=")[[1]][1])
      # the part to the left of , is the deterministic part + random bit in last term
      eqn.det <- str_trim(strsplit(eqn, ";")[[1]][1])
      # the eqn for the errors is after the ,
      eqn.stoc <- str_trim(strsplit(eqn, ";")[[1]][2])
      # the name is before the ~
      err.name <- str_trim(strsplit(eqn.stoc, "~")[[1]][1])
      # check that the name of the errors matches before and after the ,
      tmp <- strsplit(eqn.det, "[+]")[[1]] # tmp[length(tmp)] should be error part
      tmp <- strsplit(tmp[length(tmp)], "[*]")[[1]] # tmp[length(tmp)] should be err.name
      # these 2 lines get the bit between ( and ,
      eqn.stoc.mean <- strsplit(eqn.stoc, ",")[[1]][1]
      eqn.stoc.mean <- str_trim(strsplit(eqn.stoc.mean, "[(]")[[1]][2])
      # these 2 lines get the bit between , and )
      eqn.stoc.var <- strsplit(eqn.stoc, "[)]")[[1]][1]
      eqn.stoc.var <- str_trim(strsplit(eqn.stoc.var, ",")[[1]][2])

      # Set up and test the time designation for the equation term
      if (!str_detect(eqn.name, "_")) stop(paste("toLatex.marssMODEL: in the attribute, equation, ", eqn.name, " is missing the time designation.  should look like _{t}.", sep = ""))
      eqn.name.time <- str_trim(strsplit(eqn.name, "_")[[1]][2])
      eqn.name <- str_trim(strsplit(eqn.name, "_")[[1]][1])
      eqn.labs <- attr(x, paste(toupper(eqn.name), ".names", sep = ""))
      # get rid of the bit before the =
      eqn.det <- str_trim(strsplit(eqn.det, "=")[[1]][2])

      # Set up and test the time designation for error term
      if (!str_detect(err.name, "_")) stop(paste("toLatex.marssMODEL: in the attribute, equation, ", err.name, " is missing the time designation.  should look like _{t}.", sep = ""))
      err.name.time <- str_trim(strsplit(err.name, "_")[[1]][2])
      err.name <- str_trim(strsplit(err.name, "_")[[1]][1])


      #############################################################
      # Build the latex for equation
      #############################################################
      #########################################
      # start with the eqn name (left of =)
      #########################################
      if (!time.expand) {
        time.sub <- eqn.name.time
      } else {
        time.sub <- paste("{", as.character(t), "}", sep = "")
      }
      mat.body <- paste(paste(eqn.name, "_{", eqn.labs, "}", sep = ""), collapse = "\\\\\n")
      mat.tex <- paste(mat.head, "\n", mat.body, "\n", mat.foot, "_", time.sub, "=", sep = "")
      eqn.tex <- c(eqn.tex, mat.tex)
      headfoot <- list(eqn.name = eqn.name, mat.head = mat.head, mat.foot = mat.foot)

      #############################################################
      # send the deterministic part to right of = as char string like Bx_t+U
      # and build the latex for the equation
      ##############################################################
      if (time.expand) {
        eqn.tex <- c(eqn.tex, build.eqn.tex(eqn.det, eqn.special, marssMODEL, greek, digits, simplify, headfoot, t = t))
      } else {
        eqn.tex <- c(eqn.tex, build.eqn.tex(eqn.det, eqn.special, marssMODEL, greek, digits, simplify, headfoot))
      }
      # put the stochastic part on it's own line
      eqn.tex <- c(eqn.tex, "\\\\")

      ##############################################################
      # build the stochastic part of the equation
      ##############################################################
      # start with the eqn name (left of ~)
      if (!time.expand) {
        time.sub <- err.name.time
      } else {
        time.sub <- paste("{", as.character(t), "}", sep = "")
      }
      mat.body <- paste(paste(err.name, "_{", eqn.labs, "}", sep = ""), collapse = "\\\\\n")
      mat.tex <- paste(mat.head, "\n", mat.body, "\n", mat.foot, "_", time.sub, " \\sim", sep = "")
      eqn.tex <- c(eqn.tex, mat.tex)
      # add the distribution designation. will be in front of (
      tmp <- strsplit(eqn.stoc, "[(]")[[1]][1]
      tmp <- str_trim(strsplit(tmp, "~")[[1]][2])
      # wrap in pmatrix
      eqn.tex <- c(eqn.tex, tmp, "\\begin{pmatrix}")

      # add the mean
      if (!time.expand) {
        eqn.tex <- c(eqn.tex, build.eqn.tex(eqn.stoc.mean, eqn.special, marssMODEL, greek, digits, simplify, headfoot))
      } else {
        eqn.tex <- c(eqn.tex, build.eqn.tex(eqn.stoc.mean, eqn.special, marssMODEL, greek, digits, simplify, headfoot, t = t))
      }

      # separate from variance with a ,
      eqn.tex <- c(eqn.tex, ", ")
      # add the variance
      if (!time.expand) {
        eqn.tex <- c(eqn.tex, build.eqn.tex(eqn.stoc.var, eqn.special, marssMODEL, greek, digits, simplify, headfoot))
      } else {
        eqn.tex <- c(eqn.tex, build.eqn.tex(eqn.stoc.var, eqn.special, marssMODEL, greek, digits, simplify, headfoot, t = t))
      }

      # close
      eqn.tex <- c(eqn.tex, "\\end{pmatrix}")

      # end of equation i
      eqn.tex <- c(eqn.tex, "\\\\\n \\\\")
    } # end of equation (like x=..., w ~ MVN(...)
    # get rid of trailing two \\\\ after the last eqn
    eqn.tex <- eqn.tex[-1 * length(eqn.tex)]
    # add the header and footer for the equation
    all.t.eqn.tex <- c(all.t.eqn.tex, eqn.head, eqn.tex, eqn.foot)
  } # over all t for which we should print

  if(output=="rawtex") return(all.t.eqn.tex)

  ##############################################################
  # Create the latex file
  ##############################################################
  if (length(math.sty)) {
    sty <- paste("\\usepackage{", math.sty, "}", sep = "")
  }
  sty <- c(sty, paste("\\usepackage[", orientation, "]{geometry}", sep = ""))
  if ("Rnw" %in% output) {
    sty <- c(sty, "\\usepackage{Sweave}")
  }

  # Set up the file name
  if (is.null(file) & replace) tmp <- "marssMODEL"
  if (is.null(file) & !replace) tmp <- tempfile("marssMODEL", tmpdir = ".")
  if (!is.null(file)) tmp <- file
  # If the user put on a extension of .tex, .dvi or .pdf remove that
  if (tolower(str_sub(tmp, -4)) == ".tex") tmp <- str_sub(tmp, -4)
  if (tolower(str_sub(tmp, -4)) == ".pdf") tmp <- str_sub(tmp, -4)
  if (tolower(str_sub(tmp, -4)) == ".dvi") tmp <- str_sub(tmp, -4)

  tmptex.filename <- paste(tmp, "tex", sep = ".")
  cat("\\documentclass{report}", sty, "\\begin{document}\\pagestyle{empty}",
    all.t.eqn.tex, "\\end{document}\n",
    file = tmptex.filename, sep = "\n"
  )

  if ("pdf" %in% output) latexcmd <- "pdflatex"
  if ("dvi" %in% output) latexcmd <- "latex"
  if ("pdf" %in% output | "dvi" %in% output) {
    # sys(paste(latexcmd,  "-interaction=scrollmode", shQuote(tmp)), output = FALSE)
    system(paste(latexcmd, shQuote(tmp)))
    system(paste("rm", paste(str_replace(tmp, "[\\]", "/"), "aux", sep = ".")))
    system(paste("rm", paste(str_replace(tmp, "[\\]", "/"), "log", sep = ".")))
  }
  if ("Rnw" %in% output) {
    system(paste("cp", paste(str_replace(tmp, "[\\]", "/"), "tex", sep = "."), paste(str_replace(tmp, "[\\]", "/"), "Rnw", sep = ".")))
  }
  if (!("tex" %in% output)) {
    system(paste("rm", paste(str_replace(tmp, "[\\]", "/"), "tex", sep = ".")))
  }
} # end of toLatex.marssMODEL

get.mat.tex <- function(x, greek, digits) {
  # called by build.eqn.tex
  # x is a matrix for the parameter
  # use Hmisc util to convert to tabular form
  mat <- array(list(), dim = dim(x))
  tmp.fun <- function(x, digits = 2, greek = greek) {
    ifelse(is.numeric(x[[1]]), round(x[[1]], digits = digits), Hmisc::latexTranslate(x[[1]], greek = greek))
  }
  # this will be a list matrix whether or not x is a list matrix
  mat[, ] <- lapply(x, tmp.fun, digits = digits, greek = greek)
  # now deal with any +,*,-1* in the matrices
  tmp.fun.2 <- function(x) {
    if (is.character(x[[1]])) {
      x <- str_replace_all(x, "[+][-]", "-")
      x <- str_replace_all(x, "[-]1[*]", "-")
      x <- str_replace_all(x, "[+]1[*]", "+")
      x <- str_replace_all(x, "[*]", "$\\\\\\times$")
      x <- str_replace_all(x, "[+]", "$+$")
      x <- str_replace_all(x, "[-]", "$-$")
      # get rid of all the $ because the matrix will be wrapped in equation
      x <- str_replace_all(x, "[$]", " ")
    } else {
      x[[1]]
    }
  }
  # mat is still a list matrix
  mat[, ] <- lapply(mat, tmp.fun.2)

  # now use the Hmisc util to translate to tabular form
  tmp.tex <- Hmisc::latexTabular(mat)[[1]]
  # latexTabular is replacing _ with \\\\_
  tmp.tex <- str_replace_all(tmp.tex, "[\\][\\]_", "_")
  tmp.tex <- strsplit(tmp.tex, "\n")[[1]]
  tmp.tex <- tmp.tex[-1]
  tmp.tex <- tmp.tex[-1 * length(tmp.tex)]
  # tmp.tex <- tmp.tex[-1 * length(tmp.tex)]
  tmp.tex[length(tmp.tex)] <- str_replace_all(tmp.tex[length(tmp.tex)], "[\\]", "")

  tmp.body <- paste(tmp.tex, collapse = "\n")
  return(tmp.body)
}

build.eqn.tex <- function(eqn, eqn.special, x, greek, digits, simplify, headfoot, t = NULL) {
  # if simplify is FALSE, then will need t if parameters are time-varying
  # eqn.special are the terms like x, y, v, w
  # greek means whether to turn things like alpha into \alpha
  tmp.eqn.tex <- c()
  # get the model from x (x is a marssMODEL object)
  fixed <- x$fixed
  free <- x$free
  model.dims <- attr(x, "model.dims")
  m <- model.dims[["x"]][1]
  n <- model.dims[["y"]][1]
  model.el <- attr(x, "par.names")
  eqn.name <- headfoot$eqn.name
  mat.head <- headfoot$mat.head
  mat.foot <- headfoot$mat.foot
  # this takes a text form of an equation and builds the tex
  for (el in strsplit(eqn, "[+]")[[1]]) {
    show.el <- TRUE
    for (el2 in strsplit(el, "[*]")[[1]]) {
      show.el2 <- TRUE
      # if el2 has a time subscript, get the time (usually t or t-1 but might be anything)
      # el2.time is the time subscript and el2 is the parameter name
      if (str_detect(el2, "_")) {
        if (missing(t)) { # then not time.expand
          el2.time <- str_trim(strsplit(el2, "_")[[1]][2])
        } else {
          el2.time <- paste("{", as.character(t), "}", sep = "")
        }
        el2 <- strsplit(el2, "_")[[1]][1]
      } else {
        el2.time <- NULL
      }
      el2 <- str_trim(el2) # strip any extra while space
      # Detect if el2 is a number, like 0, and print as a number
      options(warn = -1) # don't print warnings if el2 cannot be converted to numeric
      if (!is.na(as.numeric(el2))) el2 <- as.numeric(el2)
      options(warn = 0)

      if (!is.numeric(el2) & !(el2 %in% c(model.el, eqn.special))) stop("toLatex.marssMODEL: in the attribute, equation, one of the names does not match a name in fixed or model.dims.")
      # if els is a number, just print as a number
      if (is.numeric(el2)) {
        mat.tex <- as.character(el2)
      } else {
        if (el2 %in% eqn.special) {
          ######################################
          # if el2 is something like x, y, v or w
          ######################################
          # if el2 is a eqn.name then use its names, otherwise it is a error name
          if (el2 %in% eqn.special[seq(1, length(eqn.special), 2)]) {
            el2.names <- el2
          } else {
            el2.names <- eqn.special[which(eqn.special == el2) - 1]
          }
          labs <- attr(x, paste(toupper(el2.names), ".names", sep = ""))
          mat.body <- paste(paste(el2, "_{", labs, "}", sep = ""), collapse = "\\\\\n")
          if (is.null(el2.time)) stop(paste("toLatex.marssMODEL: in the attribute, equation, ", el2, " is missing the time designation.  should look like _{t}.", sep = ""))
          mat.tex <- paste(mat.head, "\n", mat.body, "\n", mat.foot, "_", el2.time, sep = "")
        } else {
          ######################################
          # el2 is a parameter name
          ######################################
          if (!is.null(colnames(free[[el2]]))) { # this adds the parameter name to the element name, like B.1
            colnames(free[[el2]]) <- paste(el2, colnames(free[[el2]]), sep = ".")
          }
          if (model.dims[[el2]][3] == 1) {
            mat <- fixed.free.to.formula(fixed[[el2]], free[[el2]], model.dims[[el2]][1:2])
          } else {
            mat <- array(NA, dim = model.dims[[el2]])
            for (i in 1:model.dims[[el2]][3]) {
              mat[, , i] <- fixed.free.to.formula(fixed[[el2]], free[[el2]], model.dims[[el2]][1:2])
            }
          }
          # if simplify=TRUE only show mats that are not all 0 and that are not identity
          if (simplify) show.el <- show.el & !(all(sapply(mat, identical, 0)))
          if (simplify) {
            if (length(dim(mat)) == 3) show.el2 <- !all(apply(mat, 3, is.identity))
            if (length(dim(mat)) == 2) show.el2 <- !is.identity(mat)
          }
          if (show.el2 & show.el) {
            # if mat is 2D then it is not time-varying so no time subscript
            if (!is.matrix(mat)) { # then it is 3D
              # if mat is 3D but all the mats are equal, then it is not time-varying
              # The following tests for that
              mat2D <- array(mat, dim = c(dim(mat)[1] * dim(mat)[2], dim(mat)[3])) # reform to be 2D
              # test equality within row; will return TRUE if dim 3 = 1
              if (all(apply(mat, 1, vector.all.equal))) { # not time-varying
                mat <- sub3D(mat, t = 1) # not time-varying so just use the mat at t=1
              } else { # is time-varying
                if (!simplify) { # then there will be an equation for each parameter at time t
                  if (missing(t)) stop("build.eqn.tex: if simplify is FALSE and parameters are time-varying, need to specify t")
                  mat <- sub3D(mat, t = t)
                } else { # do simplify
                  # step 1: figure out which parameter elements are time-varying
                  mat2D <- array(mat, dim = c(dim(mat)[1] * dim(mat)[2], dim(mat)[3])) # reform to be 2D, matrix(vec(mat),t)
                  tv.elem <- !apply(mat2D, 1, vector.all.equal) # which elements are time varying
                  # step 2: replace that element with the name el2(i)_t
                  tv.name <- paste(el2, "(", 1:dim(mat2D)[1], ")_", el2.time, sep = "")
                  # step 3: set up a mat to be t=1; will replace time=varying ones
                  mat <- sub3D(mat, t = 1)
                  # replace any time-varying elements with tv.name
                  mat[tv.elem] <- tv.name[tv.elem]
                }
              }
            }
            mat.body <- get.mat.tex(mat, greek, digits)
            mat.tex <- paste(mat.head, "\n", mat.body, "\n", mat.foot, sep = "")
          }
        }
      }
      if (show.el & show.el2) tmp.eqn.tex <- c(tmp.eqn.tex, mat.tex)
    }
    if (show.el) tmp.eqn.tex <- c(tmp.eqn.tex, "+")
  }
  tmp.eqn.tex <- tmp.eqn.tex[-1 * length(tmp.eqn.tex)]
  return(tmp.eqn.tex)
}

parameters.time.varying <- function(eqn, eqn.special, x) {
  # this will search through the equation parameters and determine if any parameters are time-varying
  # get the model from x (x is a marssMODEL object)
  fixed <- x$fixed
  free <- x$free
  if (!is.null(attr(x, "model.dims"))) {
    model.dims <- attr(x, "model.dims")
  } else {
    model.dims <- x$model.dims
  }
  model.el <- names(fixed)

  # this finds the bits between + in the equation
  for (el in strsplit(eqn, "[+]")[[1]]) {
    for (el2 in strsplit(el, "[*]")[[1]]) { # this finds the bits between * in el
      if (str_detect(el2, "_")) {
        el2.time <- str_trim(strsplit(el2, "_")[[1]][2])
        el2 <- strsplit(el2, "_")[[1]][1]
      } else {
        el2.time <- NULL # no time subscript
        # the parameter should not be time-varying then....
      }
      el2 <- str_trim(el2) # strip any extra while space
      # Detect if el2 is a number, like 0, and print as a number
      options(warn = -1) # don't print warnings if el2 cannot be converted to numeric
      if (!is.na(as.numeric(el2))) {
        el2 <- as.numeric(el2)
        if (!is.null(el2.time)) stop("toLatex.marssMODEL: equation has a number with a time-subscript.  That doesn't make sense.\nCheck the equation attribute for the marssMODEL.\n")
      }
      options(warn = 0)

      if (!is.numeric(el2) & !(el2 %in% c(model.el, eqn.special))) stop("toLatex.marssMODEL: in the attribute, equation, one of the names does not match a name in fixed or model.dims.")
      # if els is a number, just print as a number
      if (is.numeric(el2)) {
        next # go to next parameter in model
      }
      if (el2 %in% eqn.special) {
        next # ignore the specials.  They always are time-varying
      }
      if (model.dims[[el2]][3] == 1) next # it is not time-varying
      # Even if 3D, it might not be time-varying if all elements are the same
      testtv <- FALSE
      if (!all(apply(fixed[[el2]], 1, vector.all.equal))) testtv <- TRUE # time-varying
      mat2D <- array(free[[el2]], dim = c(dim(free[[el2]])[1] * dim(free[[el2]])[2], dim(free[[el2]])[3]))
      if (!all(apply(mat2D, 1, vector.all.equal))) testtv <- TRUE # time-varying
      if (!testtv) { # not time-varying
        next
      } else {
        if (is.null(el2.time)) stop(paste("toLatex.marssMODEL: something is wrong. ", el2, " in the equation attributes has no time subscript, yet it is time-varying in the marssMODEL object.\nCannot make latex because the correct time-subscript label is unknown.\n", sep = ""))
        return(TRUE)
      }
    } # el2 the bits between * for el's
  } # el the bits between +
  return(FALSE)
} # end of function
