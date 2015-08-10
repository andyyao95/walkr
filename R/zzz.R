.onLoad <- function(libname, pkgname) {
  
  ## compress the vignette PDF to fix CMD Check WARNING
  
  Sys.setenv("_R_BUILD_COMPACT_VIGNETTES_" = "gs+qpdf")
}