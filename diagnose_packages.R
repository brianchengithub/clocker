#!/usr/bin/env Rscript
# ============================================================================
# 'clocker' diagnostic script
# ============================================================================
#
# Run this when something isn't working. It prints:
#   - Which optional dependencies are installed and at what version
#   - Whether each upstream backend (sesame, EpiDISH, DunedinPACE,
#     methylCIPHER, EpiMitClocks) loads its data correctly
#   - A live test of clocker on a tiny synthetic input, exercising the
#     imputation, sex-inference, and clock-projection paths
#   - The current state of clocker's caches and reference_betas resolution
#
# Usage:
#   Rscript diagnose_packages.R   > diagnostics.log 2>&1
#
# Share the resulting log when reporting issues.

cat("=== clocker diagnostics ===\n\n")
cat("R version: ", R.version.string, "\n")
cat("Platform:  ", R.version$platform, "\n")
cat("Date/Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# Local helper (clocker may not be loadable yet)
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------------------
# Helper to print package status compactly
# ---------------------------------------------------------------------------
report_pkg <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %-20s OK     v%s\n", pkg,
                as.character(packageVersion(pkg))))
    TRUE
  } else {
    cat(sprintf("  %-20s NOT INSTALLED\n", pkg))
    FALSE
  }
}


# ---------------------------------------------------------------------------
# 1. clocker itself
# ---------------------------------------------------------------------------
cat("--- clocker ---\n")
clocker_ok <- report_pkg("clocker")
if (clocker_ok) {
  cat("  Exported functions: ",
      paste(ls(asNamespace("clocker"))[!grepl("^\\.", ls(asNamespace("clocker")))][1:10],
            collapse = ", "), " ...\n")
} else {
  cat("  Install with: devtools::install_github(\"<your-username>/clocker\")\n")
}


# ---------------------------------------------------------------------------
# 2. SeSAMe (IDAT preprocessing)
# ---------------------------------------------------------------------------
cat("\n--- SeSAMe ---\n")
sesame_ok <- report_pkg("sesame")
if (sesame_ok) {
  ns <- asNamespace("sesame")
  fns <- ls(ns)
  cat("  Has openSesame: ", "openSesame" %in% fns, "\n")
  cat("  Has getBetas:   ", "getBetas"   %in% fns, "\n")
  if ("openSesame" %in% fns) {
    cat("  openSesame signature:\n    ")
    print(args(sesame::openSesame))
  }
}
report_pkg("sesameData")


# ---------------------------------------------------------------------------
# 3. EpiDISH (cell-type deconvolution)
# ---------------------------------------------------------------------------
cat("\n--- EpiDISH ---\n")
if (report_pkg("EpiDISH")) {
  cat("  Trying to load centDHSbloodDMC.m...\n")
  tryCatch({
    e <- new.env(parent = emptyenv())
    utils::data("centDHSbloodDMC.m", package = "EpiDISH", envir = e)
    if (exists("centDHSbloodDMC.m", envir = e)) {
      ref <- get("centDHSbloodDMC.m", envir = e)
      cat(sprintf("  OK: %d probes x %d cell types: %s\n",
                  nrow(ref), ncol(ref),
                  paste(colnames(ref), collapse = ", ")))
    } else {
      cat("  WARNING: dataset loaded but object not found\n")
    }
  }, error = function(e) cat("  ERROR:", e$message, "\n"))
}


# ---------------------------------------------------------------------------
# 4. DunedinPACE
# ---------------------------------------------------------------------------
cat("\n--- DunedinPACE ---\n")
if (report_pkg("DunedinPACE")) {
  ns <- asNamespace("DunedinPACE")
  cat("  Has PACEProjector: ", "PACEProjector" %in% ls(ns), "\n")
  if ("PACEProjector" %in% ls(ns)) {
    cat("  PACEProjector signature:\n    ")
    print(args(DunedinPACE::PACEProjector))
  }
}


# ---------------------------------------------------------------------------
# 5. methylCIPHER
# ---------------------------------------------------------------------------
cat("\n--- methylCIPHER ---\n")
if (report_pkg("methylCIPHER")) {
  ns <- asNamespace("methylCIPHER")
  fns <- ls(ns)
  calc <- fns[grepl("^calc", fns)]
  cat("  ", length(calc), " calc* functions available\n")
  if (length(calc)) cat("    e.g.: ", paste(head(calc, 8), collapse = ", "), "\n")

  cat("  Datasets shipped with methylCIPHER:\n")
  pkg_data <- utils::data(package = "methylCIPHER")
  if (nrow(pkg_data$results) > 0) {
    items <- pkg_data$results[, "Item"]
    cat("    ", paste(head(items, 12), collapse = ", "),
        if (length(items) > 12) sprintf(" ... (+%d)", length(items) - 12) else "",
        "\n")
  } else {
    cat("    (none — coefficient files may be missing)\n")
  }

  # Try a representative load
  for (d in c("Horvath1_CpGs", "PCClocks_CpGs", "AdaptAge_CpGs")) {
    e <- new.env(parent = emptyenv())
    res <- tryCatch({
      utils::data(list = d, package = "methylCIPHER", envir = e)
      if (exists(d, envir = e)) "OK" else "MISSING"
    }, error = function(err) paste("ERROR:", conditionMessage(err)),
       warning = function(w) paste("WARN:", conditionMessage(w)))
    cat(sprintf("    load %-20s -> %s\n", d, res))
  }
}


# ---------------------------------------------------------------------------
# 6. EpiMitClocks
# ---------------------------------------------------------------------------
cat("\n--- EpiMitClocks ---\n")
if (report_pkg("EpiMitClocks")) {
  for (d in c("dataETOC3", "estETOC3", "epiTOCcpgs3", "EpiCMITcpgs", "Replitali")) {
    e <- new.env(parent = emptyenv())
    res <- tryCatch({
      utils::data(list = d, package = "EpiMitClocks", envir = e)
      if (length(ls(e)) > 0L) sprintf("OK (%s)", paste(ls(e), collapse = ", "))
      else "MISSING"
    }, error = function(err) paste("ERROR:", conditionMessage(err)),
       warning = function(w) paste("WARN:", conditionMessage(w)))
    cat(sprintf("    load %-20s -> %s\n", d, res))
  }
}


# ---------------------------------------------------------------------------
# 7. clocker's optional helpers
# ---------------------------------------------------------------------------
cat("\n--- Optional helpers (clocker has fallbacks for all of these) ---\n")
for (p in c("matrixStats", "qs2", "digest", "progress", "curl")) report_pkg(p)


# ---------------------------------------------------------------------------
# 8. clocker cache + reference_betas resolution
# ---------------------------------------------------------------------------
cat("\n--- clocker paths and resolution ---\n")
if (clocker_ok) {
  cat("  options(clocker.cache_dir):       ",
      getOption("clocker.cache_dir") %||% "(unset)", "\n")
  cat("  Sys.getenv('CLOCKER_CACHE'):      ",
      ifelse(nzchar(Sys.getenv("CLOCKER_CACHE")), Sys.getenv("CLOCKER_CACHE"), "(unset)"),
      "\n")
  cat("  R_user_dir('clocker', 'cache'):   ",
      tools::R_user_dir("clocker", "cache"), "\n")
  cat("  options(clocker.reference_betas): ",
      getOption("clocker.reference_betas") %||% "(unset)", "\n")
  cat("  Sys.getenv('CLOCKER_REFERENCE_BETAS'): ",
      ifelse(nzchar(Sys.getenv("CLOCKER_REFERENCE_BETAS")),
             Sys.getenv("CLOCKER_REFERENCE_BETAS"), "(unset)"),
      "\n")
  ext_path <- system.file("extdata", "reference_betas.rds", package = "clocker")
  cat("  inst/extdata/reference_betas.rds: ",
      if (nzchar(ext_path)) ext_path else "(not found)", "\n")
}


# ---------------------------------------------------------------------------
# 9. Live end-to-end test on synthetic data
# ---------------------------------------------------------------------------
cat("\n--- Live end-to-end test ---\n")
if (clocker_ok) {
  # Build a tiny synthetic beta matrix with some NAs
  set.seed(42)
  n_probes <- 1000L
  n_samples <- 6L
  probes <- paste0("cg", sprintf("%08d", seq_len(n_probes)))
  betas <- matrix(rbeta(n_probes * n_samples, 2, 5),
                   nrow = n_probes, ncol = n_samples,
                   dimnames = list(probes, paste0("S", seq_len(n_samples))))
  # Inject missingness
  betas[sample(n_probes, 30), 1] <- NA
  betas[sample(n_probes, 200), 2] <- NA   # zero-shot path

  cat(sprintf("  Synthetic input: %d probes x %d samples (NAs: %d)\n",
              n_probes, n_samples, sum(is.na(betas))))
  cat("  Calling clocker(verbose = TRUE)...\n\n")

  res <- tryCatch(
    clocker::clocker(betas, verbose = TRUE),
    error = function(e) {
      cat("  ERROR: ", e$message, "\n")
      NULL
    }
  )

  if (!is.null(res)) {
    cat(sprintf("\n  -> Output: %d rows x %d columns\n", nrow(res), ncol(res)))
    cat("  -> Columns: ",
        paste(head(colnames(res), 12), collapse = ", "),
        if (ncol(res) > 12) sprintf(" ... (+%d)", ncol(res) - 12) else "", "\n")
  }
} else {
  cat("  Skipped: clocker not installed\n")
}

cat("\n=== diagnostics complete ===\n")
cat("Share this log when reporting issues.\n")
