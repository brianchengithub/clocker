# ============================================================================
# Main Entry Point: clocker()
#
# clocker() is the renamed primary entry. Aliases are preserved for back-
# compat with callers that used the older names from the prior package
# version (then called quickclocks):
#   * qclock()              -> deprecated, kept working
#   * calculate_clocks()    -> deprecated alias
#   * run_epigenetic_clocks() -> deprecated alias
# ============================================================================


#' Calculate epigenetic clocks from DNA methylation data
#'
#' Computes 40+ epigenetic clocks from raw IDAT files or pre-processed beta
#' matrices. Performs kNN imputation against your reference_betas.rds,
#' cell-type deconvolution (EpiDISH), methylQC-style sex inference, and
#' clock projection in a single call.
#'
#' @section Imputation:
#' Missing values are imputed in two passes:
#' \enumerate{
#'   \item For samples with less than \code{knn_zero_shot_threshold}
#'   missing-probe fraction, the k most similar OTHER samples in the same
#'   input batch are found by Pearson correlation, and missing probes are
#'   imputed as a Gaussian-kernel-weighted mean across those k neighbours.
#'   \item Samples above the threshold (or single-sample runs, or any
#'   probes the kNN neighbours also lack) fall through to zero-shot
#'   imputation using the per-probe means in
#'   \code{reference_betas.rds}.
#' }
#'
#' @section Sex inference:
#' Uses the exact methylQC sex caller (curated chrX/chrY probe panels,
#' threshold sweep, cluster regression, +/- 5 sigma orthogonal-distance
#' bands). When SigDFs are available (IDAT input), per-sample signals are
#' total intensity (MG + MR + UG + UR). For beta-matrix input, median
#' beta over the same probes is used as a proxy. See \code{infer_sex()}.
#'
#' @section Age fallback:
#' If \code{Age} is missing from pheno (and Age is needed by PC-Clocks),
#' Horvath2 is used as the placeholder, then Horvath1 if Horvath2 is
#' unavailable. The user is alerted whenever a fallback is used.
#'
#' @param input One of: a beta matrix (CpGs x samples), a path to an .rds /
#'   .qs2 / .csv beta-matrix file, a directory containing IDAT files, or a
#'   character vector of .idat file paths.
#' @param pheno Optional data.frame with sample-level phenotypes. Must have
#'   rownames matching colnames(betas), or a sample_id-style column.
#'   Recognized columns:
#'   \itemize{
#'     \item \code{Age} (numeric)
#'     \item \code{Female} (0/1) OR \code{Sex} ("M"/"F"/"Male"/"Female",
#'       any case). Coerced to integer 0 (Male) / 1 (Female) / NA via
#'       \code{coerce_female_indicator()}.
#'   }
#' @param n_cores Number of CPU cores. NULL = auto-detect.
#' @param reference_path Optional path to your \code{reference_betas.rds}.
#'   If NULL, resolved by \code{resolve_reference_betas_path()} in this
#'   order: \code{options("clocker.reference_betas")},
#'   \code{Sys.getenv("CLOCKER_REFERENCE_BETAS")},
#'   then \code{system.file("extdata", "reference_betas.rds",
#'   package = "clocker")}.
#' @param knn_k Number of nearest in-batch neighbours per missing-value
#'   imputation. Default 10.
#' @param knn_zero_shot_threshold Per-sample fraction of probes missing
#'   above which kNN is bypassed and reference probe means are used
#'   directly. Default 0.10 (10%).
#' @param missing_report_path Optional path to write per-sample, per-clock
#'   probe-coverage CSV (long form). A wide-format companion CSV is also
#'   written.
#' @param verbose Print progress.
#'
#' @return Data frame with one row per sample. Columns include:
#'   \describe{
#'     \item{Clock_*}{First-generation clocks (Horvath1, Horvath2, Hannum,
#'       PhenoAge, DNAmTL, Lin, Zhang, ...).}
#'     \item{PC_*}{Principal-component-based clocks.}
#'     \item{CellType_RPC_*, CellType_CP_*}{Cell composition fractions.}
#'     \item{InferredSex, chrX_signal, chrY_signal, InferredSexFlag,
#'       InferredSexScale}{Sex inference outputs.}
#'     \item{epiTOC2_TNSC, EpiCMIT, miAge, ...}{Mitotic clocks.}
#'     \item{DunedinPACE}{Pace of aging estimate.}
#'   }
#'
#' @examples
#' \dontrun{
#' # Beta matrix with no pheno: sex and age are inferred
#' results <- clocker(betas)
#'
#' # Provide pheno explicitly (Sex auto-coerced)
#' pheno <- data.frame(Age = c(45, 67, 22),
#'                     Sex = c("Female", "M", "f"),
#'                     row.names = c("S1", "S2", "S3"))
#' results <- clocker(betas, pheno = pheno)
#'
#' # Explicit reference_betas path
#' results <- clocker(betas,
#'                     reference_path = "~/mydata/reference_betas.rds")
#'
#' # Write missingness report
#' results <- clocker(betas, missing_report_path = "missingness.csv")
#'
#' # IDAT directory
#' results <- clocker("/path/to/idat_dir/")
#' }
#'
#' @seealso \code{\link{set_sex_reference}}, \code{\link{coerce_female_indicator}}
#' @export
clocker <- function(input,
                     pheno                   = NULL,
                     n_cores                 = NULL,
                     reference_path          = NULL,
                     knn_k                   = 10L,
                     knn_zero_shot_threshold = 0.10,
                     missing_report_path     = NULL,
                     verbose                 = TRUE) {

  start_time <- Sys.time()
  if (verbose) {
    message("=== clocker: Epigenetic Clock Calculator ===")
  }

  # ---- Stage 1: Load + validate input ----
  betas <- load_input_data(input, n_cores, verbose)
  betas <- validate_betas(betas)

  # ---- Stage 2: Platform detection / EPICv2 normalization ----
  platform <- detect_array_platform(rownames(betas))
  if (verbose) message("  Platform detected: ", platform)

  if (platform %in% c("EPICv2", "EPICv2/EPIC+", "MSA")) {
    betas <- normalize_epicv2_probes(betas, verbose = verbose)
  }

  # ---- Stage 3: Resource planning ----
  n_cores <- determine_optimal_cores(ncol(betas), nrow(betas), n_cores, verbose)
  if (verbose) message("  Using ", n_cores, " core(s)")

  # ---- Stage 4: kNN imputation against reference_betas.rds ----
  na_count <- sum(is.na(betas))
  if (na_count > 0L) {
    if (verbose) {
      log_msg("  Missing values: %d (%.2f%%); running kNN imputation...",
              na_count, 100 * na_count / length(betas), verbose = TRUE)
    }
    betas <- perform_knn_imputation(
      betas,
      reference_path      = reference_path,
      k                   = knn_k,
      zero_shot_threshold = knn_zero_shot_threshold,
      verbose             = verbose
    )
  } else if (verbose) {
    message("  No missing values to impute")
  }

  # ---- Stage 5: Compute clocks ----
  results <- compute_all_clocks(betas, pheno = pheno,
                                  n_cores = n_cores, verbose = verbose)

  # ---- Stage 6: Age acceleration ----
  results <- compute_age_acceleration(results, pheno, verbose)

  # ---- Stage 7: Optional missing-probe report ----
  if (!is.null(missing_report_path)) {
    write_missing_probe_report(missing_report_path, verbose = verbose)
  }

  # ---- Stage 8: Format output ----
  results <- format_output(results)

  if (verbose) {
    elapsed <- difftime(Sys.time(), start_time, units = "secs")
    message(sprintf("=== Done in %.1f sec ===", as.numeric(elapsed)))
  }

  results
}


# ---- Backward-compat aliases ----

#' @rdname clocker
#' @export
qclock <- function(...) {
  .Deprecated("clocker", package = "clocker")
  clocker(...)
}

#' @rdname clocker
#' @export
calculate_clocks <- function(...) {
  .Deprecated("clocker", package = "clocker")
  clocker(...)
}

#' @rdname clocker
#' @export
run_epigenetic_clocks <- function(...) {
  .Deprecated("clocker", package = "clocker")
  clocker(...)
}


# ---------------------------------------------------------------------------
# Age acceleration
# ---------------------------------------------------------------------------

#' @keywords internal
compute_age_acceleration <- function(results, pheno = NULL, verbose = TRUE) {

  age <- NULL
  if (!is.null(pheno) && "Age" %in% colnames(pheno)) {
    age <- as.numeric(pheno$Age[match(results$sample_id, rownames(pheno))])
    if (!all(is.na(age)) && verbose) {
      message("  Computing age acceleration vs. chronological Age...")
    }
  }
  if (is.null(age) || all(is.na(age))) return(results)

  age_clocks <- intersect(c("Horvath1", "Horvath2", "Hannum", "PhenoAge",
                             "Lin", "Zhang", "Zhang2019",
                             "PCHorvath1", "PCHorvath2", "PCHannum",
                             "PCPhenoAge", "PCGrimAge"),
                           colnames(results))

  for (cl in age_clocks) {
    pred <- as.numeric(results[[cl]])
    if (length(unique(stats::na.omit(age))) >= 2L) {
      fit <- tryCatch(stats::lm(pred ~ age, na.action = stats::na.exclude),
                       error = function(e) NULL)
      if (!is.null(fit)) {
        results[[paste0("Accel_", cl)]] <- as.numeric(stats::residuals(fit))
        next
      }
    }
    results[[paste0("Accel_", cl)]] <- pred - age
  }
  results
}


# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

#' @keywords internal
format_output <- function(results) {

  rename_map <- c(
    Horvath1   = "Clock_Horvath1",
    Horvath2   = "Clock_Horvath2",
    Hannum     = "Clock_Hannum",
    PhenoAge   = "Clock_PhenoAge",
    DNAmTL     = "Clock_DNAmTL",
    Lin        = "Clock_Lin",
    Zhang      = "Clock_Zhang",
    Zhang2019  = "Clock_Zhang2019",
    Bocklandt  = "Clock_Bocklandt",
    Weidner    = "Clock_Weidner",
    VidalBralo = "Clock_VidalBralo",
    Garagnani  = "Clock_Garagnani",
    AdaptAge   = "Clock_AdaptAge",
    CausAge    = "Clock_CausAge",
    DamAge     = "Clock_DamAge",
    SystemsAge = "Clock_SystemsAge",
    PCHorvath1 = "PC_Horvath1",
    PCHorvath2 = "PC_Horvath2",
    PCHannum   = "PC_Hannum",
    PCPhenoAge = "PC_PhenoAge",
    PCGrimAge  = "PC_GrimAge",
    PCDNAmTL   = "PC_DNAmTL",
    DunedinPACE = "DunedinPACE"
  )
  for (oldn in names(rename_map)) {
    if (oldn %in% colnames(results)) {
      colnames(results)[colnames(results) == oldn] <- rename_map[[oldn]]
    }
  }

  preferred <- c("sample_id",
                  "InferredSex", "InferredSexFlag", "InferredSexScale",
                  "chrX_signal", "chrY_signal")
  preferred <- intersect(preferred, colnames(results))
  others <- setdiff(colnames(results), preferred)
  results[, c(preferred, others), drop = FALSE]
}
