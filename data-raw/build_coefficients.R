# =============================================================================
# data-raw/build_coefficients.R
#
# Run this once to build inst data file:
#   data/clock_coefficients.rda
#
# This file embeds clock coefficients in the package so that:
#   * runtime does not depend on methylCIPHER or EpiMitClocks being installed
#   * the package can be used in offline / air-gapped environments
#   * the exact coefficient version used is pinned and reproducible
#
# Usage:
#   Rscript data-raw/build_coefficients.R
#
# Requirements (only at build-time):
#   methylCIPHER, EpiMitClocks
# =============================================================================

stopifnot(requireNamespace("methylCIPHER", quietly = TRUE))

clock_coefficients <- list()

# ---- methylCIPHER ----
mc_datasets <- c(
  "Horvath1_CpGs", "Horvath2_CpGs", "Hannum_CpGs", "PhenoAge_CpGs",
  "DNAmTL_CpGs", "Lin_CpGs", "Zhang_10_CpG", "Zhang2019_CpGs",
  "Bocklandt_CpG", "Weidner_CpGs", "VidalBralo_CpGs", "Garagnani_CpG",
  "PCClocks_CpGs", "HorvathOnlineRef",
  "AdaptAge_CpGs", "CausAge_CpGs", "DamAge_CpGs", "SystemsAge_CpGs",
  "EpiToc_CpGs", "EpiToc2_CpGs", "hypoClock_CpGs", "MiAge_CpGs"
)

for (d in mc_datasets) {
  cat("Loading methylCIPHER::", d, "...\n", sep = "")
  e <- new.env(parent = emptyenv())
  res <- tryCatch({
    utils::data(list = d, package = "methylCIPHER", envir = e)
    if (exists(d, envir = e)) {
      clock_coefficients[[d]] <- get(d, envir = e)
      "ok"
    } else "not_found"
  }, error = function(err) paste("error:", conditionMessage(err)))
  cat("  -> ", res, "\n", sep = "")
}

# ---- EpiMitClocks (optional) ----
if (requireNamespace("EpiMitClocks", quietly = TRUE)) {
  epi_datasets <- list(
    "dataETOC3"       = "dataETOC3.l",
    "estETOC3"        = "estETOC3.m",
    "epiTOCcpgs3"     = "epiTOCcpgs3.v",
    "cugpmitclockCpG" = "cugpmitclockCpG.v",
    "EpiCMITcpgs"     = "epiCMIT.df",
    "Replitali"       = c("replitali.coe", "replitali.cpg.v")
  )
  for (load_name in names(epi_datasets)) {
    cat("Loading EpiMitClocks::", load_name, "...\n", sep = "")
    e <- new.env(parent = emptyenv())
    tryCatch({
      utils::data(list = load_name, package = "EpiMitClocks", envir = e)
      for (obj in epi_datasets[[load_name]]) {
        if (exists(obj, envir = e)) {
          clock_coefficients[[obj]] <- get(obj, envir = e)
        }
      }
    }, error = function(err) cat("  -> error:", conditionMessage(err), "\n"))
  }
}

# ---- Save ----
out_path <- file.path("data", "clock_coefficients.rda")
dir.create("data", showWarnings = FALSE)
save(clock_coefficients, file = out_path, compress = "xz", compression_level = 9)
cat("\nSaved", length(clock_coefficients), "coefficient datasets to", out_path, "\n")
cat("File size: ", round(file.info(out_path)$size / 1024, 1), " KB\n", sep = "")
