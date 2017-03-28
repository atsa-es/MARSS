## Set up env for globals
.onLoad <- function(libname, pkgname) {
  assign("pkg_globals", new.env(), envir=parent.env(environment()))
  assign("alldefaults", 1, pkg_globals)
  assign("kem.methods", 1, pkg_globals)
  assign("allowed.methods", 1, pkg_globals)
  assign("common.allowed.in.MARSS.call", 1, pkg_globals)
}