#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Set default options for the package
  op <- options()
  op_quickclocks <- list(
    quickclocks.verbose = TRUE,
    quickclocks.cache_dir = file.path(path.expand("~"), ".epigenetic_clock_cache")
  )
  toset <- !(names(op_quickclocks) %in% names(op))
  if (any(toset)) options(op_quickclocks[toset])
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "quickclocks v", utils::packageVersion(pkgname),
    " - Unified Epigenetic Clock Calculator"
  )
  packageStartupMessage(
    "Run calculate_clocks() to get started."
  )
}
