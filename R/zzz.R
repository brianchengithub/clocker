# ============================================================================
# Package Load and Attach Hooks
# ============================================================================


.onLoad <- function(libname, pkgname) {
  # Set default package options
  op <- options()
  op_quickclocks <- list(
    quickclocks.verbose = TRUE,
    quickclocks.cache_dir = file.path(Sys.getenv("HOME"), ".epigenetic_clock_calculator")
  )
  toset <- !(names(op_quickclocks) %in% names(op))
  if (any(toset)) options(op_quickclocks[toset])

  invisible()
}


.onAttach <- function(libname, pkgname) {
  version <- tryCatch(
    as.character(utils::packageVersion(pkgname)),
    error = function(e) "2.0.0"
  )
  packageStartupMessage(
    sprintf("quickclocks v%s - Unified Epigenetic Clock Calculator", version),
    "\nRun calculate_clocks() to get started."
  )
}
