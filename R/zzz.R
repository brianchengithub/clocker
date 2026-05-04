# ============================================================================
# Package Load and Attach Hooks
# ============================================================================

# Private package environment - replaces .GlobalEnv writes throughout the
# package. Coefficient datasets, lazy-loaded data, and side-effecting reads
# from upstream packages are all routed here so the user's workspace is not
# polluted.
.qc_env <- new.env(parent = emptyenv())


.onLoad <- function(libname, pkgname) {
  op <- options()
  op_clocker <- list(
    clocker.verbose   = TRUE,
    clocker.cache_dir = NULL,   # overrides default cache location if set
    clocker.data_dir  = NULL    # overrides default reference DB location
  )
  toset <- !(names(op_clocker) %in% names(op))
  if (any(toset)) options(op_clocker[toset])

  invisible()
}


.onAttach <- function(libname, pkgname) {
  version <- tryCatch(
    as.character(utils::packageVersion(pkgname)),
    error = function(e) "2.1.0"
  )
  packageStartupMessage(sprintf(
    "clocker v%s -- Unified Epigenetic Clock Calculator\nRun clocker() to get started.",
    version))
}
