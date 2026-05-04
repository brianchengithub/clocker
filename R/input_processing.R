# ============================================================================
# Input Processing
#
# Handles IDAT file loading via SeSAMe and beta matrix preparation.
# When loading IDATs, also computes per-sample chrX/chrY *total signal
# intensities* from the SigDFs (the methylQC sex caller's preferred input)
# and caches them in .qc_env so SigDFs themselves can be released after
# beta extraction. Beta-only inputs skip this step; sex inference then
# falls back to a beta-proxy that uses the same downstream algorithm.
# ============================================================================


#' Process IDAT files using SeSAMe
#'
#' Pipeline applied: "QCDPB" (QC, Detection p-value, Pval drop, Background,
#' BMIQ). SigDFs are processed streaming so we never keep the full SigDF
#' list in memory:
#'   1. Run openSesame to get one SigDF per sample
#'   2. For each SigDF: compute sex intensities (chrX, chrY) into a small
#'      data.frame, cache via cache_sex_signals()
#'   3. Extract betas via getBetas(sdf) and cbind into a matrix
#'   4. Drop SigDFs
#'
#' Single-sample handling (FIX H6): when a single .idat path is supplied,
#' the underlying sample basename is preserved in colnames(betas) rather
#' than being replaced by an arbitrary index.
#'
#' @keywords internal
process_idat_files <- function(input, n_cores, verbose) {

  if (!requireNamespace("sesame", quietly = TRUE)) {
    stop("sesame package required for IDAT processing. ",
         "Install with: BiocManager::install('sesame')")
  }

  if (length(input) == 1L && dir.exists(input)) {
    if (verbose) message("  Processing IDAT directory: ", input)
    sdfs <- sesame::openSesame(input, prep = "QCDPB", func = NULL, BPPARAM = NULL)
  } else {
    if (verbose) message("  Processing ", length(input), " IDAT file(s)...")
    paths <- input
    if (length(paths) == 1L) {
      base <- sub("_(Grn|Red)\\.idat(\\.gz)?$", "", paths, ignore.case = TRUE)
      sample_name <- basename(base)
      sdfs <- sesame::openSesame(base, prep = "QCDPB", func = NULL)
      if (is.list(sdfs)) names(sdfs) <- sample_name
    } else {
      sdfs <- sesame::openSesame(paths, prep = "QCDPB", func = NULL)
    }
  }

  # Wrap single-SigDF returns in a list for uniform handling
  if (!is.list(sdfs) || (!is.null(sdfs$Probe_ID) && !is.list(sdfs[[1]]))) {
    sdfs <- list(sdfs)
    names(sdfs) <- "Sample_1"
  }
  if (is.null(names(sdfs))) {
    names(sdfs) <- paste0("Sample_", seq_along(sdfs))
  }

  # ---- Compute sex intensities from SigDFs (methylQC preferred input) ----
  sex_intensities <- tryCatch(
    compute_sex_signals_from_sdfs(sdfs),
    error = function(e) {
      if (verbose) message("    Sex intensity calc failed: ", e$message)
      NULL
    })
  if (!is.null(sex_intensities)) {
    cache_sex_signals(sex_intensities, scale = "intensity")
    if (verbose) {
      n_ok <- sum(!is.na(sex_intensities$chrX) & !is.na(sex_intensities$chrY))
      log_msg("  Cached methylQC-style intensity signals (%d/%d samples)",
              n_ok, nrow(sex_intensities), verbose = verbose)
    }
  }

  # ---- Extract betas ----
  betas <- do.call(cbind, lapply(sdfs, sesame::getBetas))

  if (is.null(colnames(betas)) || any(colnames(betas) == "")) {
    colnames(betas) <- names(sdfs)
  }

  # SigDFs no longer needed -- free memory
  rm(sdfs)
  invisible(gc(verbose = FALSE))

  betas
}


#' Load beta values from various input formats
#' @keywords internal
load_input_data <- function(input, n_cores, verbose) {

  # Reset cached sex signals from any previous run
  .qc_env$sex_signals <- NULL

  if (is.matrix(input) || is.data.frame(input)) {
    betas <- as.matrix(input)
  } else if (is.character(input) && length(input) >= 1L) {
    is_idat <- grepl("\\.idat(\\.gz)?$", input, ignore.case = TRUE) |
               (length(input) == 1L && dir.exists(input))
    if (any(is_idat)) {
      betas <- process_idat_files(input, n_cores, verbose)
    } else if (length(input) == 1L) {
      f <- input
      if (grepl("\\.qs2?$", f) && requireNamespace("qs2", quietly = TRUE)) {
        betas <- qs2::qs_read(f)
      } else if (grepl("\\.rds$", f, ignore.case = TRUE)) {
        betas <- readRDS(f)
      } else if (grepl("\\.csv$", f, ignore.case = TRUE)) {
        betas <- as.matrix(utils::read.csv(f, row.names = 1, check.names = FALSE))
      } else if (grepl("\\.t(ab|sv)$", f, ignore.case = TRUE)) {
        betas <- as.matrix(utils::read.table(f, header = TRUE, sep = "\t",
                                                row.names = 1, check.names = FALSE))
      } else {
        stop("Unsupported file format: ", f)
      }
    } else {
      stop("Multiple file paths only supported for .idat input")
    }
  } else {
    stop("Unsupported input type: ", class(input)[1])
  }

  if (!is.matrix(betas)) betas <- as.matrix(betas)
  betas
}


#' Apply kNN imputation (and stash diagnostic frame for the orchestrator)
#' @keywords internal
perform_knn_imputation <- function(betas,
                                     reference_path,
                                     k,
                                     zero_shot_threshold,
                                     verbose) {

  imp <- knn_impute(
    betas,
    reference_path      = reference_path,
    k                   = k,
    zero_shot_threshold = zero_shot_threshold,
    verbose             = verbose
  )
  .qc_env$imputation_info <- imp$sample_info
  imp$betas
}
