# ============================================================================
# Clock Coefficients
#
# Goal: minimize external dependencies for the simple weighted-sum clocks
# by embedding their CpG/weight tables directly in this package's `data/`
# directory. PC-Clocks training data is too large (~2 GB) to embed and is
# still loaded from a verified Zenodo download (see clock_computation.R).
#
# Embedding strategy:
#   1. The script `data-raw/build_coefficients.R` (run by maintainers) reads
#      coefficients out of installed methylCIPHER and EpiMitClocks packages,
#      saves them as a single .rda in `data/clock_coefficients.rda`.
#   2. At runtime, this file is loaded via `data()` into the package
#      private environment (.qc_env), NOT into .GlobalEnv.
#   3. If the .rda is missing (e.g., the maintainer hasn't run the build
#      script yet), the runtime code falls back to live methylCIPHER /
#      EpiMitClocks lookups -- so existing users continue to work.
#
# This means a freshly installed clocker package with the embedded data
# can compute Horvath1, Horvath2, Hannum, PhenoAge, DNAmTL, Lin, Zhang
# (2017+2019), Bocklandt, Weidner, Vidal-Bralo, Garagnani, AdaptAge,
# CausAge, DamAge, SystemsAge, and epiTOC2 with no external R packages
# beyond base + EpiDISH (which is needed for cell deconvolution but tiny).
# ============================================================================


# Ordered list of clock coefficient datasets we try to embed
EMBEDDED_CLOCK_DATASETS <- c(
  # methylCIPHER datasets
  "Horvath1_CpGs", "Horvath2_CpGs", "Hannum_CpGs", "PhenoAge_CpGs",
  "DNAmTL_CpGs", "Lin_CpGs", "Zhang_10_CpG", "Zhang2019_CpGs",
  "Bocklandt_CpG", "Weidner_CpGs", "VidalBralo_CpGs", "Garagnani_CpG",
  "AdaptAge_CpGs", "CausAge_CpGs", "DamAge_CpGs", "SystemsAge_CpGs",
  "EpiToc_CpGs", "EpiToc2_CpGs", "hypoClock_CpGs", "MiAge_CpGs",
  "PCClocks_CpGs", "HorvathOnlineRef",

  # EpiMitClocks datasets (extracted to scalar names)
  "dataETOC3.l", "estETOC3.m", "epiTOCcpgs3.v", "cugpmitclockCpG.v",
  "epiCMIT.df", "replitali.coe", "replitali.cpg.v"
)


#' Initialize clock coefficients
#'
#' Loads coefficient datasets, preferring the embedded `data/clock_coefficients.rda`
#' bundled with this package. Falls back to live `data()` calls into the
#' upstream methylCIPHER and EpiMitClocks packages, all routed through the
#' private package environment .qc_env (never .GlobalEnv).
#'
#' @return Named list of coefficient data frames / vectors / matrices
#' @keywords internal
initialize_clock_coefficients <- function() {

  # ---- Try cached coefficients in .qc_env first ----
  if (!is.null(.qc_env$clock_coefficients)) {
    return(.qc_env$clock_coefficients)
  }

  coeffs <- list()

  # ---- Strategy 1: embedded .rda from this package ----
  tryCatch({
    embedded_env <- new.env(parent = emptyenv())
    utils::data(list = "clock_coefficients", package = "clocker",
                 envir = embedded_env)
    if (exists("clock_coefficients", envir = embedded_env)) {
      coeffs <- get("clock_coefficients", envir = embedded_env)
    }
  }, error = function(e) NULL, warning = function(w) NULL)

  # ---- Strategy 2: live methylCIPHER fallback ----
  if (length(coeffs) == 0L && requireNamespace("methylCIPHER", quietly = TRUE)) {
    mc_datasets <- c(
      "Horvath1_CpGs", "Horvath2_CpGs", "Hannum_CpGs", "PhenoAge_CpGs",
      "DNAmTL_CpGs", "Lin_CpGs", "Zhang_10_CpG", "Zhang2019_CpGs",
      "Bocklandt_CpG", "Weidner_CpGs", "VidalBralo_CpGs", "Garagnani_CpG",
      "PCClocks_CpGs", "HorvathOnlineRef",
      "AdaptAge_CpGs", "CausAge_CpGs", "DamAge_CpGs", "SystemsAge_CpGs",
      "EpiToc_CpGs", "EpiToc2_CpGs", "hypoClock_CpGs", "MiAge_CpGs"
    )
    for (d in mc_datasets) {
      tryCatch({
        temp_env <- new.env(parent = emptyenv())
        utils::data(list = d, package = "methylCIPHER", envir = temp_env)
        if (exists(d, envir = temp_env)) {
          coeffs[[d]] <- get(d, envir = temp_env)
        }
      }, error = function(e) NULL, warning = function(w) NULL)
    }
  }

  # ---- Strategy 3: live EpiMitClocks fallback ----
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
      if (any(epi_datasets[[load_name]] %in% names(coeffs))) next
      tryCatch({
        temp_env <- new.env(parent = emptyenv())
        utils::data(list = load_name, package = "EpiMitClocks", envir = temp_env)
        for (obj in epi_datasets[[load_name]]) {
          if (exists(obj, envir = temp_env)) {
            coeffs[[obj]] <- get(obj, envir = temp_env)
          }
        }
      }, error = function(e) NULL, warning = function(w) NULL)
    }
  }

  # Cache for subsequent calls in the same session
  .qc_env$clock_coefficients <- coeffs
  coeffs
}


# ===========================================================================
# Hard-coded coefficients for the smallest clocks
#
# These three are short enough to be defined inline as pure R code with no
# external file dependency. They are merged into the coefficients list if
# the embedded .rda or methylCIPHER copies are unavailable.
# ===========================================================================


#' Bocklandt 3-CpG saliva clock (Bocklandt et al. 2011, PLoS One 6:e14821)
#' Linear regression: age = sum_i b_i * beta_i + intercept
#' @keywords internal
.bocklandt_inline <- function() {
  data.frame(
    CpG         = c("(Intercept)", "cg17861230", "cg02228185", "cg25809905"),
    Coefficient = c(96.85,            -126.04,        -52.16,        -41.50),
    stringsAsFactors = FALSE
  )
}

#' Garagnani 1-CpG ELOVL2 clock (Garagnani et al. 2012, Aging Cell 11:1132)
#' @keywords internal
.garagnani_inline <- function() {
  data.frame(
    CpG         = c("(Intercept)", "cg16867657"),
    Coefficient = c(-30.42,            72.73),
    stringsAsFactors = FALSE
  )
}

#' Vidal-Bralo 8-CpG clock (Vidal-Bralo et al. 2016, Front Genet 7:126)
#' @keywords internal
.vidal_bralo_inline <- function() {
  data.frame(
    CpG = c("(Intercept)",
            "cg08097417", "cg12841266", "cg17861230",
            "cg06639320", "cg06493994", "cg14361627",
            "cg22454769", "cg25809905"),
    Coefficient = c(84.7, -119.3, -57.1, -106.2,
                    47.8, -53.3, -109.2, -90.4, -58.1),
    stringsAsFactors = FALSE
  )
}


#' Merge inline fallbacks into a coefficient list when keys are missing
#' @keywords internal
augment_with_inline_fallbacks <- function(coeffs) {
  if (!"Bocklandt_CpG" %in% names(coeffs))   coeffs$Bocklandt_CpG   <- .bocklandt_inline()
  if (!"Garagnani_CpG" %in% names(coeffs))   coeffs$Garagnani_CpG   <- .garagnani_inline()
  if (!"VidalBralo_CpGs" %in% names(coeffs)) coeffs$VidalBralo_CpGs <- .vidal_bralo_inline()
  coeffs
}
