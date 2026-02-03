# ============================================================================
# Input Processing Functions
# IDAT loading, reference beta loading, imputation
# ============================================================================


#' Load IDAT files from directory using SeSAMe
#' @param idat_dir Path to directory containing IDAT files
#' @param verbose Print progress
#' @return Beta matrix (probes as rows, samples as columns)
#' @keywords internal
load_idat_directory <- function(idat_dir, verbose = TRUE) {

  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("Package 'sesame' is required for IDAT processing.\n",
         "Install with: BiocManager::install('sesame')")
  }

  # Find all IDAT files
  idat_files <- list.files(idat_dir, pattern = "_Grn\\.idat$|_Red\\.idat$",
                           recursive = TRUE, full.names = TRUE, ignore.case = TRUE)

  if (length(idat_files) == 0) {
    stop("No IDAT files found in: ", idat_dir)
  }

  # Get unique sample prefixes
  sample_prefixes <- unique(gsub("_(Grn|Red)\\.idat$", "", idat_files, ignore.case = TRUE))
  log_msg("Found %d samples", length(sample_prefixes), verbose = verbose)

  # Process with SeSAMe openSesame
  log_msg("Processing with SeSAMe (NOOB normalization, dye-bias correction, pOOBAH)...",
          verbose = verbose)

  betas <- sesame::openSesame(
    idat_dir,
    prep = "QCDPB",
    func = sesame::getBetas
  )

  # openSesame returns a named vector for 1 sample, matrix for multiple
  # Convert to matrix if needed (probes as rows, samples as columns)
  if (is.numeric(betas) && !is.matrix(betas)) {
    probe_names <- names(betas)
    sample_name <- if (length(sample_prefixes) == 1) {
      basename(sample_prefixes[1])
    } else {
      "Sample1"
    }
    betas <- matrix(betas, ncol = 1, dimnames = list(probe_names, sample_name))
  }

  # If openSesame returned samples as rows (samples < probes), transpose
  if (!is.null(betas) && is.matrix(betas) && nrow(betas) < ncol(betas)) {
    betas <- t(betas)
  }

  return(betas)
}


#' Load reference betas for imputation
#'
#' Searches multiple locations for the reference betas file used during
#' imputation. Falls back gracefully to row-median imputation if not found.
#'
#' @param verbose Print progress
#' @return Named numeric vector of reference beta values, or NULL if not found
#' @keywords internal
load_reference_betas <- function(verbose = TRUE) {

  # Helper to safely check if a file exists
  safe_file_exists <- function(p) {
    tryCatch({
      if (is.null(p) || !is.character(p) || length(p) != 1) return(FALSE)
      if (is.na(p) || nchar(p) == 0) return(FALSE)
      file.exists(p)
    }, error = function(e) FALSE)
  }

  # Build list of candidate paths
  ref_paths <- character(0)

  # 1. Installed package location
  tryCatch({
    p <- system.file("extdata", "reference_betas.rds", package = "quickclocks")
    if (nchar(p) > 0) ref_paths <- c(ref_paths, p)
  }, error = function(e) NULL)

  # 2. Development: relative to current working directory
  ref_paths <- c(ref_paths, file.path(getwd(), "inst", "extdata", "reference_betas.rds"))

  # 3. Try relative to source file location (only works when sourced)
  tryCatch({
    src_file <- sys.frame(1)$ofile
    if (!is.null(src_file) && is.character(src_file) && nchar(src_file) > 0) {
      ref_paths <- c(ref_paths,
        file.path(dirname(src_file), "..", "inst", "extdata", "reference_betas.rds"),
        file.path(dirname(src_file), "..", "..", "inst", "extdata", "reference_betas.rds")
      )
    }
  }, error = function(e) NULL)

  # 4. Try to find via installed package path
  tryCatch({
    pkg_path <- find.package("quickclocks", quiet = TRUE)
    if (length(pkg_path) > 0 && nchar(pkg_path[1]) > 0) {
      ref_paths <- c(ref_paths, file.path(pkg_path[1], "extdata", "reference_betas.rds"))
    }
  }, error = function(e) NULL)

  # 5. Also check for the original filename
  ref_paths <- c(ref_paths,
    file.path(getwd(), "inst", "extdata", "final_cg_means_01032026.rds"),
    file.path(getwd(), "final_cg_means_01032026.rds")
  )

  # Check each path safely
  for (path in ref_paths) {
    if (safe_file_exists(path)) {
      tryCatch({
        ref <- readRDS(path)
        log_msg("Loaded %d reference probe values", length(ref), verbose = verbose)
        return(ref)
      }, error = function(e) NULL)
    }
  }

  log_msg("Reference betas file not found. Imputation will use row medians only.",
          verbose = verbose)
  return(NULL)
}


#' Perform smart imputation of missing beta values
#'
#' Uses reference betas when available, falls back to row medians,
#' then to 0.5 for any remaining NAs.
#'
#' @param betas Beta matrix (probes as rows, samples as columns)
#' @param reference_betas Named vector of reference beta values (or NULL)
#' @param verbose Print progress
#' @return Imputed beta matrix
#' @keywords internal
perform_smart_imputation <- function(betas, reference_betas, verbose = TRUE) {

  n_missing_before <- sum(is.na(betas))

  if (n_missing_before == 0) {
    log_msg("No missing values to impute", verbose = verbose)
    return(betas)
  }

  log_msg("Missing values: %d (%.2f%%)",
          n_missing_before,
          100 * n_missing_before / length(betas),
          verbose = verbose)

  if (!is.null(reference_betas)) {
    probes_with_missing <- rownames(betas)[apply(betas, 1, function(x) any(is.na(x)))]

    for (probe in probes_with_missing) {
      na_idx <- is.na(betas[probe, ])
      if (probe %in% names(reference_betas)) {
        betas[probe, na_idx] <- reference_betas[probe]
      } else {
        row_median <- median(betas[probe, ], na.rm = TRUE)
        betas[probe, na_idx] <- if (!is.na(row_median)) row_median else 0.5
      }
    }
  } else {
    for (i in seq_len(nrow(betas))) {
      na_idx <- is.na(betas[i, ])
      if (any(na_idx)) {
        row_median <- median(betas[i, ], na.rm = TRUE)
        betas[i, na_idx] <- if (!is.na(row_median)) row_median else 0.5
      }
    }
  }

  n_missing_after <- sum(is.na(betas))
  log_msg("Imputed %d values", n_missing_before - n_missing_after, verbose = verbose)

  return(betas)
}
