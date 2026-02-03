# ============================================================================
# Main Entry Point
# The single exported function: calculate_clocks()
# ============================================================================


#' Calculate Epigenetic Clocks
#'
#' Main function to calculate 40+ epigenetic clocks from DNA methylation data.
#' Automatically detects whether input is a directory path (IDAT files) or
#' a numeric matrix of beta values.
#'
#' @param input Either:
#'   \itemize{
#'     \item A character string path to directory containing IDAT files (can be nested)
#'     \item A path to an RDS or CSV file containing a beta matrix
#'     \item A numeric matrix of beta values with CpG probe IDs as rownames
#'           and sample IDs as colnames
#'   }
#' @param pheno Optional data.frame with sample phenotype data. Should have:
#'   \itemize{
#'     \item Row names matching sample IDs in beta matrix
#'     \item "Age" column (numeric, chronological age in years)
#'     \item "Female" column (numeric: 1=female, 0=male, 0.5=unknown)
#'   }
#'   Required for accurate PC clock calculations.
#' @param n_cores Number of CPU cores to use. Default (NULL) automatically selects
#'   optimal number based on input size, available cores, and RAM.
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return A data.frame with sample IDs as rows and clock values as columns.
#'   Includes cell composition estimates, sex inference, and all available clocks.
#'
#' @examples
#' \dontrun{
#' # From IDAT directory
#' results <- calculate_clocks("/path/to/idat/directory")
#'
#' # From beta matrix
#' results <- calculate_clocks(my_beta_matrix)
#'
#' # With phenotype data for PC clocks
#' pheno <- data.frame(Age = c(45, 52, 38), Female = c(1, 0, 1),
#'                     row.names = colnames(my_beta_matrix))
#' results <- calculate_clocks(my_beta_matrix, pheno = pheno)
#' }
#'
#' @export
calculate_clocks <- function(input, pheno = NULL, n_cores = NULL, verbose = TRUE) {

  start_time <- Sys.time()

  # ============================================================
  # Print banner
  # ============================================================
  if (verbose) {
    cat("\n")
    cat("==============================================================\n")
    cat("       UNIFIED EPIGENETIC CLOCK CALCULATOR v2.0.0\n")
    cat("==============================================================\n")
    cat("\n")
    cat("Integrating clocks from:\n")
    cat("  - SeSAMe (IDAT preprocessing)\n")
    cat("  - EpiDISH (cell type deconvolution: RPC + CP methods)\n")
    cat("  - DunedinPACE (pace of aging)\n")
    cat("  - PC-Clocks (Levine Lab)\n")
    cat("  - epiTOC2 (mitotic clocks)\n")
    cat("  - methylCIPHER (40+ clocks)\n")
    cat("\n")
  }

  # ============================================================
  # Detect input type and load data
  # ============================================================

  if (is.character(input) && length(input) == 1) {
    if (!dir.exists(input) && !file.exists(input)) {
      stop("Input path does not exist: ", input)
    }

    if (dir.exists(input)) {
      log_msg("Loading IDAT files from: %s", input, verbose = verbose)
      betas <- load_idat_directory(input, verbose = verbose)
    } else if (grepl("\\.(rds|RDS)$", input)) {
      log_msg("Loading beta matrix from RDS: %s", input, verbose = verbose)
      betas <- readRDS(input)
    } else if (grepl("\\.(csv|CSV)$", input)) {
      log_msg("Loading beta matrix from CSV: %s", input, verbose = verbose)
      betas <- as.matrix(read.csv(input, row.names = 1, check.names = FALSE))
    } else {
      stop("Unrecognized file type. Provide a directory, .rds, or .csv file.")
    }

  } else if (is.matrix(input) || is.data.frame(input)) {
    log_msg("Using provided beta matrix", verbose = verbose)
    betas <- as.matrix(input)

  } else {
    stop("Input must be either:\n",
         "  - A path to IDAT directory or beta matrix file (character string)\n",
         "  - A numeric matrix with CpG probes as rownames and sample IDs as colnames")
  }

  # Validate beta matrix
  betas <- validate_betas(betas)

  n_samples <- ncol(betas)
  n_probes <- nrow(betas)

  log_msg("Input: %d samples x %d probes", n_samples, n_probes, verbose = verbose)

  # ============================================================
  # Detect platform
  # ============================================================

  platform <- detect_array_platform(rownames(betas))
  log_msg("Detected platform: %s", platform, verbose = verbose)

  # ============================================================
  # Determine optimal thread count
  # ============================================================

  n_cores <- determine_optimal_cores(n_samples, n_probes, n_cores, verbose)
  log_msg("Using %d CPU core(s)", n_cores, verbose = verbose)

  # ============================================================
  # Load reference betas and perform imputation
  # ============================================================

  reference_betas <- tryCatch(
    load_reference_betas(verbose),
    error = function(e) {
      if (verbose) message("  Reference betas loading failed: ", e$message)
      NULL
    }
  )

  log_msg("\n--- Imputation ---", verbose = verbose)
  betas <- perform_smart_imputation(betas, reference_betas, verbose)

  # ============================================================
  # Check clock availability
  # ============================================================

  log_msg("\n--- Checking Clock Availability ---", verbose = verbose)
  availability <- check_clock_availability(rownames(betas), verbose)

  # ============================================================
  # Calculate all clocks
  # ============================================================

  log_msg("\n--- Calculating Clocks ---", verbose = verbose)
  results <- compute_all_clocks(betas, pheno, n_cores, verbose)

  # ============================================================
  # Compile final output
  # ============================================================

  if (!"sample_id" %in% colnames(results)) {
    results <- cbind(sample_id = colnames(betas), results)
  }

  rownames(results) <- results$sample_id
  results <- as.data.frame(results, stringsAsFactors = FALSE)

  # ============================================================
  # Summary
  # ============================================================

  runtime <- difftime(Sys.time(), start_time, units = "mins")

  if (verbose) {
    n_clocks <- ncol(results) - 1
    cat("\n")
    cat("==============================================================\n")
    cat("                        COMPLETE\n")
    cat("==============================================================\n")
    cat(sprintf("  Samples:  %d\n", n_samples))
    cat(sprintf("  Platform: %s\n", platform))
    cat(sprintf("  Clocks:   %d\n", n_clocks))
    cat(sprintf("  Runtime:  %.2f minutes\n", as.numeric(runtime)))
    cat("==============================================================\n")
    cat("\n")
  }

  return(results)
}


#' Wrapper for backward compatibility
#' @rdname calculate_clocks
#' @export
run_epigenetic_clocks <- calculate_clocks
