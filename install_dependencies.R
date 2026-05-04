#!/usr/bin/env Rscript
# ============================================================================
# Install Dependencies for the 'clocker' Epigenetic Clock Calculator
# ============================================================================
#
# Run once before installing/using clocker:
#
#   source("install_dependencies.R")
#
# Then install clocker itself:
#
#   devtools::install_github("<your-username>/clocker")
#   library(clocker)
#   results <- clocker("/path/to/idats")
#
# This script installs a minimal, curated set of dependencies. Optional
# packages (qs2, progress, digest, matrixStats, curl) are also installed
# because clocker uses them when available and falls back gracefully when
# they are not.

cat("============================================================\n")
cat("clocker - Dependency Installation\n")
cat("============================================================\n\n")


# ---- Helper -----------------------------------------------------------------

install_if_missing <- function(pkg, source = "CRAN", repo = NULL) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  SKIP    %-30s (already installed)\n", pkg))
    return(TRUE)
  }
  cat(sprintf("  INSTALL %-30s from %s\n", pkg, source))
  ok <- tryCatch({
    if (source == "CRAN") {
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    } else if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager",
                         repos = "https://cloud.r-project.org", quiet = TRUE)
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
    } else if (source == "GitHub") {
      if (!requireNamespace("devtools", quietly = TRUE)) {
        install.packages("devtools",
                         repos = "https://cloud.r-project.org", quiet = TRUE)
      }
      devtools::install_github(repo, quiet = TRUE, upgrade = "never")
    }
    TRUE
  }, error = function(e) {
    cat(sprintf("    WARNING: %s failed: %s\n", pkg, e$message))
    FALSE
  })
  if (ok && requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("    OK      %s\n", pkg))
    TRUE
  } else {
    cat(sprintf("    FAILED  %s\n", pkg))
    FALSE
  }
}


# ---- CRAN -------------------------------------------------------------------
# Only what clocker actually uses (or uses optionally with a fallback).
cat("\n--- CRAN packages ---\n")

cran_packages <- c(
  "devtools",       # for install_github (used to install clocker itself)
  "remotes",        # alternative to devtools
  # Optional but recommended (clocker has graceful fallbacks for all of these):
  "matrixStats",    # vectorized rowMedians / colMedians
  "qs2",            # fast binary I/O for cached manifests + PC-Clocks data
  "digest",         # SHA-256 verification of cached PC-Clocks data
  "progress",       # nicer progress bars
  "curl"            # robust downloads with retry
)

cran_results <- vapply(cran_packages,
                       function(pkg) install_if_missing(pkg, "CRAN"),
                       logical(1))


# ---- Bioconductor -----------------------------------------------------------
cat("\n--- Bioconductor packages ---\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager",
                   repos = "https://cloud.r-project.org", quiet = TRUE)
}

bioc_packages <- c(
  # Required only if you read IDAT files directly:
  "sesame",       # IDAT preprocessing (NOOB, dye bias, pOOBAH)
  "sesameData",   # SeSAMe annotation cache

  # Required for cell-type deconvolution (recommended; small):
  "EpiDISH"       # RPC + CP cell composition
)

bioc_results <- vapply(bioc_packages,
                       function(pkg) install_if_missing(pkg, "Bioconductor"),
                       logical(1))


# ---- GitHub-only packages ---------------------------------------------------
cat("\n--- GitHub packages ---\n")
cat("  (skip any you do not need — clocker degrades gracefully)\n")

github_packages <- list(
  list(pkg = "DunedinPACE",   repo = "danbelsky/DunedinPACE"),
  list(pkg = "EpiMitClocks",  repo = "aet21/EpiMitClocks"),
  list(pkg = "methylCIPHER",  repo = "MorganLevineLab/methylCIPHER")
)

github_results <- vapply(github_packages,
                         function(item)
                           install_if_missing(item$pkg, "GitHub", item$repo),
                         logical(1))


# ---- Summary ---------------------------------------------------------------
cat("\n============================================================\n")
cat("Installation summary\n")
cat("============================================================\n\n")

print_block <- function(label, names_vec, ok_vec) {
  cat(label, ":\n", sep = "")
  for (i in seq_along(names_vec)) {
    cat(sprintf("  [%-3s] %s\n",
                if (ok_vec[i]) "OK" else "FAIL", names_vec[i]))
  }
  cat("\n")
}
print_block("CRAN",          cran_packages,                              cran_results)
print_block("Bioconductor",  bioc_packages,                              bioc_results)
print_block("GitHub",        vapply(github_packages, `[[`, "", "pkg"),   github_results)


# ---- Optional: pre-cache SeSAMe annotation ---------------------------------
if (requireNamespace("sesameData", quietly = TRUE)) {
  cat("--- Pre-caching SeSAMe annotation data ---\n")
  for (plat in c("EPIC", "HM450")) {
    res <- tryCatch({
      sesameData::sesameDataCache(plat)
      TRUE
    }, error = function(e) {
      cat(sprintf("  WARNING: failed to cache %s: %s\n", plat, e$message))
      FALSE
    })
    if (isTRUE(res)) cat(sprintf("  OK  %s annotation cached\n", plat))
  }
} else {
  cat("--- SeSAMe annotation pre-caching skipped (sesameData unavailable) ---\n")
}


# ---- Final messages ---------------------------------------------------------
cat("\n============================================================\n")
all_required <- all(c(cran_results[c("devtools", "remotes")],
                       bioc_results[c("sesame", "EpiDISH")]))
if (all_required) {
  cat("Required dependencies installed. You can now run:\n")
  cat("  devtools::install_github(\"<your-username>/clocker\")\n")
  cat("  library(clocker)\n")
  cat("  results <- clocker(\"/path/to/idats\")\n")
} else {
  cat("WARNING: some required dependencies failed.\n")
  cat("Re-check the log above. Common issues:\n")
  cat("  - GitHub rate limit: set GITHUB_PAT environment variable\n")
  cat("  - Bioconductor version mismatch: BiocManager::install(version = '...')\n")
  cat("  - System libraries missing: see each package's INSTALL notes\n")
}
cat("\nNotes on optional packages:\n")
cat("  - DunedinPACE / methylCIPHER / EpiMitClocks: clocker computes only\n")
cat("    the clocks whose backing packages are available.\n")
cat("  - PC-Clocks training data (~2 GB) downloads on first use to\n")
cat("    tools::R_user_dir(\"clocker\", \"cache\") (override via\n")
cat("    options(clocker.cache_dir = ...) or env CLOCKER_CACHE).\n")
cat("============================================================\n")
