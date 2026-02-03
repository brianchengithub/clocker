#!/usr/bin/env Rscript
#' ============================================================================
#' Unified Epigenetic Clock Calculator
#' ============================================================================
#' 
#' Simple R interface for calculating 40+ epigenetic clocks from DNA methylation data.
#' 
#' @author Brian Chen
#' @version 2.0.0
#' ============================================================================


# ============================================================================
# Manifest Download and Caching Functions
# ============================================================================

#' Get manifest cache directory
#' @return Path to cache directory
#' @keywords internal
get_manifest_cache_dir <- function() {

  cache_dir <- file.path(Sys.getenv("HOME"), ".epigenetic_clock_calculator", "manifests")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  return(cache_dir)
}

#' Get manifest URL for a platform
#' @param platform Platform name (EPIC, EPICv2, HM450, etc.)
#' @return URL to manifest file
#' @keywords internal
get_manifest_url <- function(platform) {
  base_url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno"
  urls <- list(
    "MSA" = paste0(base_url, "/MSA/MSA.hg38.manifest.tsv.gz"),
    "EPICv2" = paste0(base_url, "/EPICv2/EPICv2.hg38.manifest.tsv.gz"),
    "EPIC+" = paste0(base_url, "/EPIC+/EPIC+.hg38.manifest.tsv.gz"),
    "EPIC" = paste0(base_url, "/EPIC/EPIC.hg38.manifest.tsv.gz"),
    "HM450" = paste0(base_url, "/HM450/HM450.hg38.manifest.tsv.gz")
  )
  return(urls[[platform]])
}

#' Download and cache manifest for a platform
#' @param platform Platform name
#' @param verbose Print progress
#' @return Data frame with manifest data, or NULL on failure
#' @keywords internal
download_manifest <- function(platform, verbose = TRUE) {
  cache_dir <- get_manifest_cache_dir()
  
  # Prefer qs2 format for faster loading, fall back to rds
  use_qs2 <- requireNamespace("qs2", quietly = TRUE)
  cache_file_qs2 <- file.path(cache_dir, paste0(platform, ".hg38.manifest.qs2"))
  cache_file_rds <- file.path(cache_dir, paste0(platform, ".hg38.manifest.rds"))
  
  # Check for cached file (prefer qs2, then rds)
  if (use_qs2 && file.exists(cache_file_qs2)) {
    if (verbose) message("    Loading cached manifest for ", platform, " (qs2)")
    return(qs2::qs_read(cache_file_qs2))
  } else if (file.exists(cache_file_rds)) {
    if (verbose) message("    Loading cached manifest for ", platform, " (rds)")
    manifest <- readRDS(cache_file_rds)
    # Upgrade to qs2 format if available
    if (use_qs2) {
      tryCatch({
        qs2::qs_save(manifest, cache_file_qs2)
        unlink(cache_file_rds)  # Remove old rds file
        if (verbose) message("    Upgraded cache to qs2 format")
      }, error = function(e) NULL)
    }
    return(manifest)
  }
  
  # Download
  url <- get_manifest_url(platform)
  if (is.null(url)) {
    if (verbose) message("    Unknown platform: ", platform)
    return(NULL)
  }
  
  if (verbose) message("    Downloading manifest for ", platform, "...")
  
  temp_file <- tempfile(fileext = ".tsv.gz")
  
  tryCatch({
    download.file(url, temp_file, mode = "wb", quiet = !verbose)
    
    # Read the manifest
    manifest <- read.table(gzfile(temp_file), header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE, quote = "", 
                          comment.char = "", fill = TRUE)
    
    # Cache it (prefer qs2 for ~10x faster loading)
    if (use_qs2) {
      qs2::qs_save(manifest, cache_file_qs2)
      if (verbose) message("    Cached manifest: ", nrow(manifest), " probes (qs2 format)")
    } else {
      saveRDS(manifest, cache_file_rds)
      if (verbose) message("    Cached manifest: ", nrow(manifest), " probes (rds format)")
    }
    
    unlink(temp_file)
    return(manifest)
    
  }, error = function(e) {
    if (verbose) message("    Failed to download manifest: ", e$message)
    unlink(temp_file)
    return(NULL)
  })
}

#' Get chromosome information for probes from manifest
#' @param platform Platform name
#' @param verbose Print progress
#' @return Named vector: probe_id -> chromosome
#' @keywords internal
get_probe_chromosomes <- function(platform, verbose = TRUE) {
  manifest <- download_manifest(platform, verbose)
  
  if (is.null(manifest)) return(NULL)
  
  # Find probe ID column (usually first column or "Probe_ID")
  probe_col <- NULL
  for (col in c("Probe_ID", "probeID", "IlmnID", "Name")) {
    if (col %in% colnames(manifest)) {
      probe_col <- col
      break
    }
  }
  if (is.null(probe_col)) probe_col <- colnames(manifest)[1]
  
  # Find chromosome column
  chr_col <- NULL
  for (col in c("CpG_chrm", "CHR", "chr", "Chromosome", "seqnames")) {
    if (col %in% colnames(manifest)) {
      chr_col <- col
      break
    }
  }
  
  if (is.null(chr_col)) {
    if (verbose) message("    Could not find chromosome column in manifest")
    return(NULL)
  }
  
  # Create named vector
  probe_chr <- manifest[[chr_col]]
  names(probe_chr) <- manifest[[probe_col]]
  
  return(probe_chr)
}


#' Infer sex from DNA methylation data using X/Y chromosome probes
#' @param betas Beta value matrix (CpGs as rows, samples as columns)
#' @param platform Platform name (EPIC, HM450, etc.)
#' @param verbose Print progress
#' @return Named vector of sex predictions (1 = Female, 0 = Male, 0.5 = Unknown)
#' @keywords internal
infer_sex_from_betas <- function(betas, platform = "EPIC", verbose = TRUE) {
  
  n_samples <- ncol(betas)
  sex_pred <- rep(0.5, n_samples)  # Default to unknown
  names(sex_pred) <- colnames(betas)
  
  # Get chromosome information
  probe_chr <- get_probe_chromosomes(platform, verbose = FALSE)
  
  if (is.null(probe_chr)) {
    if (verbose) message("    Could not load chromosome annotations")
    return(sex_pred)
  }
  
  # Match probes to our data
  common_probes <- intersect(names(probe_chr), rownames(betas))
  
  if (length(common_probes) == 0) {
    if (verbose) message("    No probes matched manifest")
    return(sex_pred)
  }
  
  # Get X and Y chromosome probes
  x_probes <- common_probes[probe_chr[common_probes] %in% c("chrX", "X")]
  y_probes <- common_probes[probe_chr[common_probes] %in% c("chrY", "Y")]
  
  if (verbose) {
    message("    Found ", length(x_probes), " chrX probes, ", length(y_probes), " chrY probes")
  }
  
  if (length(x_probes) < 100 || length(y_probes) < 10) {
    if (verbose) message("    Insufficient sex chromosome probes for inference")
    return(sex_pred)
  }
  
  # Calculate median beta values for X and Y chromosomes
  x_betas <- betas[x_probes, , drop = FALSE]
  y_betas <- betas[y_probes, , drop = FALSE]
  
  # Median methylation on X and Y
  x_median <- apply(x_betas, 2, median, na.rm = TRUE)
  y_median <- apply(y_betas, 2, median, na.rm = TRUE)
  
  # Sex inference logic:
  # Females (XX): Higher X methylation due to X-inactivation, very low Y signal
  # Males (XY): Lower X methylation, higher Y signal
  
  # Y chromosome median is key - males have methylation, females have noise
  # Threshold: Y median > 0.1 suggests male
  
  for (i in seq_len(n_samples)) {
    if (y_median[i] > 0.2) {
      # Clear Y signal -> Male
      sex_pred[i] <- 0
    } else if (y_median[i] < 0.1 && x_median[i] > 0.3) {
      # Low Y, higher X -> Female
      sex_pred[i] <- 1
    } else {
      # Ambiguous
      sex_pred[i] <- 0.5
    }
  }
  
  if (verbose) {
    n_male <- sum(sex_pred == 0)
    n_female <- sum(sex_pred == 1)
    n_unknown <- sum(sex_pred == 0.5)
    message("    Sex inference: ", n_female, " Female, ", n_male, " Male, ", n_unknown, " Unknown")
  }
  
  return(sex_pred)
}


# ============================================================================
# Direct Clock Implementation Functions
# These bypass problematic lazy data loading in packages
# ============================================================================

#' Initialize clock coefficients from installed packages
#' @return List of coefficient data frames
#' @keywords internal
initialize_clock_coefficients <- function() {
  coeffs <- list()
  
  # ===== methylCIPHER coefficients =====
  if (requireNamespace("methylCIPHER", quietly = TRUE)) {
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
        temp_env <- new.env()
        data(list = d, package = "methylCIPHER", envir = temp_env)
        if (exists(d, envir = temp_env)) {
          coeffs[[d]] <- get(d, envir = temp_env)
        }
      }, error = function(e) NULL, warning = function(w) NULL)
    }
  }
  
  # ===== EpiMitClocks coefficients =====
  if (requireNamespace("EpiMitClocks", quietly = TRUE)) {
    epi_datasets <- list(
      "dataETOC3" = "dataETOC3.l",
      "estETOC3" = "estETOC3.m", 
      "epiTOCcpgs3" = "epiTOCcpgs3.v",
      "cugpmitclockCpG" = "cugpmitclockCpG.v",
      "EpiCMITcpgs" = "epiCMIT.df",
      "Replitali" = c("replitali.coe", "replitali.cpg.v")
    )
    
    for (load_name in names(epi_datasets)) {
      tryCatch({
        temp_env <- new.env()
        data(list = load_name, package = "EpiMitClocks", envir = temp_env)
        obj_names <- epi_datasets[[load_name]]
        for (obj in obj_names) {
          if (exists(obj, envir = temp_env)) {
            coeffs[[obj]] <- get(obj, envir = temp_env)
          }
        }
      }, error = function(e) NULL, warning = function(w) NULL)
    }
  }
  
  return(coeffs)
}


#' Horvath age transformation (anti-log transformation)
#' Converts raw clock values to DNAm age in years
#' @param x Raw clock value (weighted sum)
#' @param adult_age Adult age constant (default 20 for Horvath clocks)
#' @return DNAm age in years
#' @keywords internal
horvath_age_transform <- function(x, adult_age = 20) {
  # Horvath uses a piecewise transformation:
  # - For ages < 20: logarithmic (raw values are log-transformed ages)
  # - For ages >= 20: linear
  # The inverse transformation:
  ifelse(x < 0,
         (adult_age + 1) * exp(x) - 1,
         (adult_age + 1) * x + adult_age)
}


#' Generic weighted sum clock calculator
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coef_df Data frame with CpG and weight columns
#' @param transform_func Optional transformation function to apply to results
#' @return Named vector of clock values
#' @keywords internal
calc_weighted_sum_clock <- function(betas, coef_df, transform_func = NULL) {
  
  # Auto-detect CpG column
  cpg_candidates <- c("CpG", "CpGmarker", "probe", "Probe", "cpg", "ID")
  cpg_col <- intersect(cpg_candidates, colnames(coef_df))[1]
  if (is.na(cpg_col)) cpg_col <- colnames(coef_df)[1]
  
  # Auto-detect weight column
  weight_candidates <- c("Coefficient", "CoefficientTraining", "weight", "Weight", "coef", "beta")
  weight_col <- intersect(weight_candidates, colnames(coef_df))[1]
  if (is.na(weight_col)) {
    # Find first numeric column that's not the CpG column
    numeric_cols <- sapply(coef_df, is.numeric)
    numeric_cols[cpg_col] <- FALSE
    weight_col <- names(which(numeric_cols))[1]
  }
  
  if (is.na(weight_col)) return(NULL)
  
  cpgs <- as.character(coef_df[[cpg_col]])
  weights <- as.numeric(coef_df[[weight_col]])
  
  # Remove intercept row if present
  intercept_mask <- cpgs %in% c("(Intercept)", "Intercept") | is.na(cpgs) | cpgs == ""
  if (any(intercept_mask)) {
    intercept <- weights[intercept_mask][1]
    if (is.na(intercept)) intercept <- 0
    cpgs <- cpgs[!intercept_mask]
    weights <- weights[!intercept_mask]
  } else {
    intercept <- 0
  }
  
  # Match CpGs
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  
  if (sum(valid) == 0) return(NULL)
  
  betas_subset <- betas[matched_idx[valid], , drop = FALSE]
  weights_valid <- weights[valid]
  
  clock_values <- intercept + colSums(betas_subset * weights_valid, na.rm = TRUE)
  
  # Apply transformation if provided (e.g., Horvath age transformation)
  if (!is.null(transform_func)) {
    clock_values <- transform_func(clock_values)
  }
  
  return(clock_values)
}


#' Calculate epiTOC2 mitotic clock directly
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named vector of TNSC values
#' @keywords internal
calc_epitoc2_direct <- function(betas, coeffs) {
  
  estETOC2 <- NULL
  
  # Try EpiMitClocks data format
  if ("dataETOC3.l" %in% names(coeffs)) {
    data_list <- coeffs[["dataETOC3.l"]]
    if (is.list(data_list) && length(data_list) >= 1) {
      estETOC2 <- data_list[[1]]
    }
  }
  
  # Try methylCIPHER data format
  if (is.null(estETOC2) && "EpiToc2_CpGs" %in% names(coeffs)) {
    estETOC2 <- coeffs$EpiToc2_CpGs
  }
  
  if (is.null(estETOC2)) return(NULL)
  
  cpgs <- rownames(estETOC2)
  if (is.null(cpgs)) return(NULL)
  
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  
  n_valid <- sum(valid)
  if (n_valid == 0) return(NULL)
  
  message(sprintf("    epiTOC2: Using %d of %d CpGs", n_valid, length(cpgs)))
  
  betas_matched <- betas[matched_idx[valid], , drop = FALSE]
  params_matched <- estETOC2[valid, , drop = FALSE]
  
  if (ncol(params_matched) >= 2) {
    delta <- params_matched[, 1]
    beta0 <- params_matched[, 2]
    
    tnsc <- apply(betas_matched, 2, function(b) {
      denom <- delta * (1 - beta0)
      denom[denom == 0] <- NA
      scores <- (b - beta0) / denom
      2 * mean(scores, na.rm = TRUE)
    })
    
    return(tnsc)
  }
  
  return(NULL)
}


#' Calculate PC Clocks directly
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named list of PC clock values
#' @keywords internal
calc_pcclocks_direct <- function(betas, coeffs) {
  
  if (!"PCClocks_CpGs" %in% names(coeffs)) return(NULL)
  
  pc_data <- coeffs$PCClocks_CpGs
  
  # PC clocks use a different methodology - need reference PCs and projections

  # This is complex and requires the full PC training data
  # For now, return NULL - we'll use the methylCIPHER function with proper data loading
  
  return(NULL)
}


#' Calculate Epigenetic Clocks
#' 
#' Main function to calculate epigenetic clocks from DNA methylation data.
#' Automatically detects whether input is a directory path (IDAT files) or 
#' a numeric matrix of beta values.
#' 
#' @param input Either:
#'   \itemize{
#'     \item A character string path to directory containing IDAT files (can be nested)
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
    # Input is a path string
    if (!dir.exists(input) && !file.exists(input)) {
      stop("Input path does not exist: ", input)
    }
    
    if (dir.exists(input)) {
      # Directory of IDAT files
      log_msg("Loading IDAT files from: %s", input, verbose = verbose)
      betas <- load_idat_directory(input, verbose = verbose)
    } else if (grepl("\\.(rds|RDS)$", input)) {
      # RDS file
      log_msg("Loading beta matrix from RDS: %s", input, verbose = verbose)
      betas <- readRDS(input)
    } else if (grepl("\\.(csv|CSV)$", input)) {
      # CSV file
      log_msg("Loading beta matrix from CSV: %s", input, verbose = verbose)
      betas <- as.matrix(read.csv(input, row.names = 1, check.names = FALSE))
    } else {
      stop("Unrecognized file type. Provide a directory, .rds, or .csv file.")
    }
    
  } else if (is.matrix(input) || is.data.frame(input)) {
    # Input is a matrix/data.frame
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
  # Load reference betas (packaged within)
  # ============================================================
  
  reference_betas <- load_reference_betas(verbose)
  
  # ============================================================
  # Perform imputation
  # ============================================================
  
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
  
  # Ensure sample_id is the first column, then convert to data.frame
  if (!"sample_id" %in% colnames(results)) {
    results <- cbind(sample_id = colnames(betas), results)
  }
  
  # Set rownames to sample IDs
  rownames(results) <- results$sample_id
  
  # Convert to data.frame (not tibble)
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  
  # ============================================================
  # Summary
  # ============================================================
  
  runtime <- difftime(Sys.time(), start_time, units = "mins")
  
  if (verbose) {
    n_clocks <- ncol(results) - 1  # Subtract sample_id column
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


# ============================================================================
# Internal Helper Functions
# ============================================================================

#' Log message with formatting
#' @keywords internal
log_msg <- function(fmt, ..., verbose = TRUE) {
  if (verbose) {
    message(sprintf(fmt, ...))
  }
}


#' Validate beta matrix
#' @keywords internal
validate_betas <- function(betas) {
  
  if (!is.numeric(betas)) {
    stop("Beta matrix must contain numeric values")
  }
  
  if (is.null(rownames(betas))) {
    stop("Beta matrix must have CpG probe IDs as rownames (e.g., 'cg00000029')")
  }
  
  if (is.null(colnames(betas))) {
    stop("Beta matrix must have sample IDs as colnames")
  }
  
  # Check that rownames look like CpG probes
  probe_pattern <- sum(grepl("^cg|^ch", rownames(betas), ignore.case = TRUE))
  if (probe_pattern < nrow(betas) * 0.5) {
    warning("Less than 50% of rownames appear to be CpG probe IDs. ",
            "Expected format: 'cg00000029', 'ch.1.1234', etc.")
  }
  
  # Check value range
  if (min(betas, na.rm = TRUE) < -0.1 || max(betas, na.rm = TRUE) > 1.1) {
    warning("Beta values outside expected range [0, 1]. ",
            "Are these M-values? Converting to beta values.")
    # Convert M-values to beta
    betas <- 2^betas / (2^betas + 1)
  }
  
  return(betas)
}


#' Detect array platform from probe IDs
#' @keywords internal
detect_array_platform <- function(probe_ids) {
  
  n_probes <- length(probe_ids)
  
  # Platform-specific probe counts (approximate)
  # MSA (Methylation Screening Array): ~285,000 probes
  # EPIC+ (EPICv2): ~930,000 probes  
  # EPIC (EPICv1): ~865,000 probes
  # 450K: ~485,000 probes
  # 27K: ~27,000 probes
  
  # Check for platform-specific probes
  has_epic_v2_probes <- any(grepl("^cg.*_TC", probe_ids)) || 
                        any(grepl("^nv_", probe_ids, ignore.case = TRUE))
  has_msa_probes <- any(grepl("^MSA", probe_ids, ignore.case = TRUE))
  
  if (has_msa_probes || (n_probes > 250000 && n_probes < 350000)) {
    return("MSA")
  } else if (has_epic_v2_probes || n_probes > 900000) {
    return("EPICv2/EPIC+")
  } else if (n_probes > 800000) {
    return("EPIC")
  } else if (n_probes > 400000) {
    return("450K")
  } else if (n_probes > 20000) {
    return("27K")
  } else {
    return("Unknown (subset)")
  }
}


#' Determine optimal number of cores based on data size and available resources
#' @keywords internal
determine_optimal_cores <- function(n_samples, n_probes, requested_cores, verbose) {
  
  # Get system info
  available_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(available_cores)) available_cores <- 1
  
  # Estimate memory requirements
  # Each sample-probe value is 8 bytes (double)
  # Plus overhead for intermediate calculations (~3x)
  bytes_per_sample <- n_probes * 8 * 3
  
  # Try to get available RAM (system-dependent)
  available_ram_gb <- tryCatch({
    if (.Platform$OS.type == "unix") {
      # Linux/Mac
      if (file.exists("/proc/meminfo")) {
        # Linux
        meminfo <- readLines("/proc/meminfo", n = 3)
        mem_free <- as.numeric(gsub("[^0-9]", "", meminfo[grep("MemAvailable|MemFree", meminfo)[1]]))
        mem_free / 1024 / 1024  # Convert KB to GB
      } else {
        # Mac - use sysctl
        mem_str <- system("sysctl -n hw.memsize", intern = TRUE)
        as.numeric(mem_str) / 1024^3  # Bytes to GB
      }
    } else {
      # Windows
      mem_str <- system("wmic OS get FreePhysicalMemory", intern = TRUE)[2]
      as.numeric(trimws(mem_str)) / 1024 / 1024  # KB to GB
    }
  }, error = function(e) {
    8  # Default assumption: 8GB available
  })
  
  if (is.na(available_ram_gb) || available_ram_gb <= 0) {
    available_ram_gb <- 8
  }
  
  # Leave 2GB for system
  usable_ram_gb <- max(available_ram_gb - 2, 1)
  
  # Calculate max samples that can be processed per core
  gb_per_sample <- bytes_per_sample / 1024^3
  max_parallel_samples <- floor(usable_ram_gb / gb_per_sample)
  
  # Optimal cores based on memory
  memory_limited_cores <- max(1, min(available_cores, floor(max_parallel_samples)))
  
  # Optimal cores based on sample count (no point using more cores than samples)
  sample_limited_cores <- min(available_cores, n_samples)
  
  # Take the minimum of all constraints
  optimal_cores <- min(memory_limited_cores, sample_limited_cores, available_cores - 1)
  optimal_cores <- max(1, optimal_cores)  # At least 1 core
  
  # If user specified cores, use that (but warn if too high)
  if (!is.null(requested_cores)) {
    if (requested_cores > optimal_cores) {
      if (verbose) {
        warning(sprintf(
          "Requested %d cores may exceed available resources. Recommended: %d cores",
          requested_cores, optimal_cores
        ))
      }
    }
    return(max(1, min(requested_cores, available_cores)))
  }
  
  return(optimal_cores)
}


#' Load reference betas from package
#' @keywords internal
load_reference_betas <- function(verbose = TRUE) {
  
  # Try multiple locations to find the reference file
  # This handles both installed package and direct sourcing scenarios
  
  ref_paths <- c(
    # 1. Installed package location (try both package names)
    system.file("extdata", "reference_betas.rds", package = "quickclocks"),
    system.file("extdata", "reference_betas.rds", package = "EpigeneticClockCalculator"),
    
    # 2. Development: relative to current working directory
    file.path(getwd(), "inst", "extdata", "reference_betas.rds")
  )
  
  # 3. Try relative to source file location (only works when sourced)
  tryCatch({
    src_file <- sys.frame(1)$ofile
    if (!is.null(src_file) && nchar(src_file) > 0) {
      ref_paths <- c(ref_paths,
        file.path(dirname(src_file), "..", "inst", "extdata", "reference_betas.rds"),
        file.path(dirname(src_file), "..", "..", "inst", "extdata", "reference_betas.rds")
      )
    }
  }, error = function(e) NULL)
  
  # 4. Try to find via installed package path
  tryCatch({
    pkg_path <- find.package("quickclocks", quiet = TRUE)
    if (length(pkg_path) > 0) {
      ref_paths <- c(ref_paths,
        file.path(pkg_path, "extdata", "reference_betas.rds")
      )
    }
  }, error = function(e) NULL)
  
  # Also check for the original filename
  ref_paths <- c(ref_paths,
    file.path(getwd(), "inst", "extdata", "final_cg_means_01032026.rds"),
    file.path(getwd(), "final_cg_means_01032026.rds")
  )
  
  # Filter out empty/NULL paths, then check each
  for (path in ref_paths) {
    tryCatch({
      if (!is.null(path) && is.character(path) && nchar(path) > 0 && file.exists(path)) {
        ref <- readRDS(path)
        log_msg("Loaded %d reference probe values", length(ref), verbose = verbose)
        return(ref)
      }
    }, error = function(e) NULL)
  }
  
  # If not found, warn user
log_msg("WARNING: Reference betas file not found.", verbose = verbose)
  log_msg("  Expected location: inst/extdata/reference_betas.rds", verbose = verbose)
  log_msg("  Imputation will use row medians only.", verbose = verbose)
  return(NULL)
}


#' Null coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Load IDAT files from directory
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
    prep = "QCDPB",  # Quality mask, Channel switch, Dye-bias, P-value mask, Background
    func = sesame::getBetas
  )
  
  # openSesame returns a named vector for 1 sample, matrix for multiple
  # Convert to matrix if needed (probes as rows, samples as columns)
  if (is.numeric(betas) && !is.matrix(betas)) {
    # Single sample: named vector -> single-column matrix
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


#' Perform smart imputation
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
  
  # Impute using reference betas if available
  if (!is.null(reference_betas)) {
    # For each probe with missing values, use reference median
    probes_with_missing <- rownames(betas)[apply(betas, 1, function(x) any(is.na(x)))]
    
    for (probe in probes_with_missing) {
      if (probe %in% names(reference_betas)) {
        na_idx <- is.na(betas[probe, ])
        betas[probe, na_idx] <- reference_betas[probe]
      } else {
        # Use row median if reference not available
        na_idx <- is.na(betas[probe, ])
        row_median <- median(betas[probe, ], na.rm = TRUE)
        if (!is.na(row_median)) {
          betas[probe, na_idx] <- row_median
        } else {
          betas[probe, na_idx] <- 0.5  # Fallback
        }
      }
    }
  } else {
    # Simple row median imputation
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


#' Check which clocks can be computed
#' @keywords internal
check_clock_availability <- function(available_probes, verbose = TRUE) {
  
  # Check for packages - some have multiple possible names
  # Note: methylCIPHER includes PC clocks
  availability <- list(
    sesame = requireNamespace("sesame", quietly = TRUE),
    dunedin_pace = requireNamespace("DunedinPACE", quietly = TRUE),
    # PC-Clocks: included in methylCIPHER
    pc_clocks = requireNamespace("methylCIPHER", quietly = TRUE) ||
                requireNamespace("methylcipher", quietly = TRUE),
    epitoc2 = requireNamespace("EpiMitClocks", quietly = TRUE) ||
              requireNamespace("epiTOC2", quietly = TRUE),
    methylcipher = requireNamespace("methylCIPHER", quietly = TRUE) ||
                   requireNamespace("methylcipher", quietly = TRUE)
  )
  
  if (verbose) {
    for (pkg in names(availability)) {
      status <- if (availability[[pkg]]) "available" else "not installed"
      # Special note for PC-Clocks
      if (pkg == "pc_clocks" && availability[[pkg]]) {
        status <- "available (via methylCIPHER)"
      }
      log_msg("  %s: %s", pkg, status, verbose = TRUE)
    }
  }
  
  return(availability)
}


#' Compute all clocks
#' @keywords internal
compute_all_clocks <- function(betas, pheno = NULL, n_cores, verbose = TRUE) {
  
  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)
  
  # ===== Direct implementations (most reliable) =====
  # These bypass the problematic lazy data loading in packages
  if (verbose) message("  Computing clocks via direct implementations...")
  
  tryCatch({
    # Initialize coefficients from packages
    coeffs <- initialize_clock_coefficients()
    
    if (length(coeffs) > 0 && verbose) {
      message("    Loaded ", length(coeffs), " coefficient datasets")
    }
    
    # ===== Direct methylCIPHER clock calculations =====
    # These are simple weighted sums
    
    mc_clocks <- list(
      "Horvath1" = "Horvath1_CpGs",
      "Horvath2" = "Horvath2_CpGs", 
      "Hannum" = "Hannum_CpGs",
      "PhenoAge" = "PhenoAge_CpGs",
      "DNAmTL" = "DNAmTL_CpGs",
      "Lin" = "Lin_CpGs",
      "Zhang" = "Zhang_10_CpG",
      "Zhang2019" = "Zhang2019_CpGs",
      "Bocklandt" = "Bocklandt_CpG",
      "Weidner" = "Weidner_CpGs",
      "VidalBralo" = "VidalBralo_CpGs",
      "Garagnani" = "Garagnani_CpG"
    )
    
    # Clocks that need Horvath age transformation (anti.trafo)
    horvath_transform_clocks <- c("Horvath1", "Horvath2")
    
    computed_direct <- c()
    
    for (clock_name in names(mc_clocks)) {
      coef_name <- mc_clocks[[clock_name]]
      if (coef_name %in% names(coeffs)) {
        tryCatch({
          # Apply Horvath age transformation for Horvath1 and Horvath2
          if (clock_name %in% horvath_transform_clocks) {
            result <- calc_weighted_sum_clock(betas, coeffs[[coef_name]], 
                                              transform_func = horvath_age_transform)
          } else {
            result <- calc_weighted_sum_clock(betas, coeffs[[coef_name]])
          }
          if (!is.null(result) && length(result) == ncol(betas)) {
            results[[clock_name]] <- as.numeric(result)
            computed_direct <- c(computed_direct, clock_name)
          }
        }, error = function(e) {
          if (verbose) message("    ", clock_name, " error: ", e$message)
        })
      }
    }
    
    # ===== epiTOC2 direct calculation =====
    tryCatch({
      toc2 <- calc_epitoc2_direct(betas, coeffs)
      if (!is.null(toc2) && length(toc2) == ncol(betas)) {
        results$epiTOC2_TNSC <- toc2
        computed_direct <- c(computed_direct, "epiTOC2_TNSC")
      }
    }, error = function(e) {
      if (verbose) message("    epiTOC2 direct error: ", e$message)
    })
    
    if (verbose && length(computed_direct) > 0) {
      message("    Direct clocks: ", paste(computed_direct, collapse = ", "))
    }
    
  }, error = function(e) {
    if (verbose) message("    Direct implementations error: ", e$message)
  })
  
  # ===== EpiDISH Cell Type Deconvolution =====
  # Use EpiDISH for cell deconvolution (best accuracy per Nature Communications 2024 benchmarking)
  # Provides both RPC (Robust Partial Correlations) and CP (Constrained Projection) methods
  if (requireNamespace("EpiDISH", quietly = TRUE)) {
    tryCatch({
      log_msg("  Estimating cell composition with EpiDISH...", verbose = verbose)
      
      # Load reference matrix - try available references in order of preference
      ref_matrix <- NULL
      ref_name <- NULL
      
      # List of references to try (in order of preference)
      # centDHSbloodDMC.m is the most commonly available (7 cell types)
      refs_to_try <- c("centDHSbloodDMC.m")
      
      for (ref in refs_to_try) {
        tryCatch({
          data(list = ref, package = "EpiDISH", envir = environment())
          ref_matrix <- get(ref)
          ref_name <- ref
          break
        }, error = function(e) NULL)
      }
      
      if (!is.null(ref_matrix)) {
        if (verbose) message("    Loaded reference: ", ref_name, " (", ncol(ref_matrix), " cell types)")
        if (verbose) message("    Cell types: ", paste(colnames(ref_matrix), collapse = ", "))
        
        # EpiDISH expects CpGs as rows, samples as columns (our format)
        # But verify and transpose if needed
        betas_for_epidish <- betas
        
        # Find overlapping probes
        common_probes <- intersect(rownames(betas_for_epidish), rownames(ref_matrix))
        if (verbose) message("    Overlapping probes with reference: ", length(common_probes))
        
        if (length(common_probes) < 100) {
          if (verbose) message("    Warning: Very few overlapping probes, cell estimates may be unreliable")
        }
        
        if (length(common_probes) > 0) {
          # Subset to common probes
          betas_subset <- betas_for_epidish[common_probes, , drop = FALSE]
          ref_subset <- ref_matrix[common_probes, , drop = FALSE]
          
          # Run RPC method (recommended - most robust to outliers)
          rpc_result <- tryCatch({
            EpiDISH::epidish(beta.m = betas_subset, ref.m = ref_subset, method = "RPC")
          }, error = function(e) {
            if (verbose) message("    EpiDISH RPC error: ", e$message)
            NULL
          })
          
          if (!is.null(rpc_result) && !is.null(rpc_result$estF)) {
            rpc_fracs <- as.data.frame(rpc_result$estF)
            
            # Ensure proper row alignment with samples
            if (nrow(rpc_fracs) == ncol(betas)) {
              # Prefix column names with method
              colnames(rpc_fracs) <- paste0("CellType_RPC_", colnames(rpc_fracs))
              
              for (col in colnames(rpc_fracs)) {
                results[[col]] <- rpc_fracs[[col]]
              }
              
              if (verbose) {
                cell_types <- gsub("CellType_RPC_", "", colnames(rpc_fracs))
                message("    EpiDISH RPC: ", length(cell_types), " cell types estimated")
              }
            }
          }
          
          # Run CP method (constrained projection - traditional Houseman)
          cp_result <- tryCatch({
            EpiDISH::epidish(beta.m = betas_subset, ref.m = ref_subset, method = "CP")
          }, error = function(e) {
            if (verbose) message("    EpiDISH CP error: ", e$message)
            NULL
          })
          
          if (!is.null(cp_result) && !is.null(cp_result$estF)) {
            cp_fracs <- as.data.frame(cp_result$estF)
            
            if (nrow(cp_fracs) == ncol(betas)) {
              # Prefix column names with method
              colnames(cp_fracs) <- paste0("CellType_CP_", colnames(cp_fracs))
              
              for (col in colnames(cp_fracs)) {
                results[[col]] <- cp_fracs[[col]]
              }
              
              if (verbose) {
                message("    EpiDISH CP: ", ncol(cp_fracs), " cell types estimated")
              }
            }
          }
        } else {
          if (verbose) message("    No overlapping probes with reference - skipping cell deconvolution")
        }
      } else {
        if (verbose) message("    Could not load any EpiDISH reference matrix")
      }
    }, error = function(e) {
      if (verbose) message("    EpiDISH deconvolution error: ", e$message)
    })
  } else {
    if (verbose) message("  EpiDISH not installed - skipping cell deconvolution")
    if (verbose) message("  Install with: BiocManager::install('EpiDISH')")
  }
  
  # ===== Sex Inference =====
  # Use custom chromosome-based sex inference function
  tryCatch({
    log_msg("  Inferring sex...", verbose = verbose)
    
    # Detect platform
    n_probes <- nrow(betas)
    platform <- "EPIC"  # Default
    if (n_probes > 900000) {
      platform <- "EPICv2"
    } else if (n_probes > 800000) {
      platform <- "EPIC"
    } else if (n_probes > 400000) {
      platform <- "HM450"
    }
    
    # Use our custom sex inference function
    sex_numeric <- infer_sex_from_betas(betas, platform = platform, verbose = verbose)
    
    if (!is.null(sex_numeric) && length(sex_numeric) == ncol(betas)) {
      # Store numeric sex (1=Female, 0=Male, 0.5=Unknown)
      results$InferredSex_Numeric <- sex_numeric
      
      # Also store character version
      sex_char <- ifelse(sex_numeric == 1, "F", ifelse(sex_numeric == 0, "M", "U"))
      results$InferredSex <- sex_char
    }
  }, error = function(e) {
    if (verbose) message("    Sex inference failed: ", e$message)
  })
  
  # ===== DunedinPACE =====
  if (requireNamespace("DunedinPACE", quietly = TRUE)) {
    tryCatch({
      log_msg("  Calculating DunedinPACE...", verbose = verbose)
      
      # DunedinPACE expects probes as rows, samples as columns (our format)
      pace <- DunedinPACE::PACEProjector(betas)
      
      if (!is.null(pace)) {
        if (is.list(pace) && "DunedinPACE" %in% names(pace)) {
          results$DunedinPACE <- pace$DunedinPACE
        } else if (is.data.frame(pace) && "DunedinPACE" %in% colnames(pace)) {
          results$DunedinPACE <- pace$DunedinPACE
        } else if (is.numeric(pace) && length(pace) == ncol(betas)) {
          results$DunedinPACE <- pace
        }
      }
    }, error = function(e) {
      if (verbose) message("    DunedinPACE failed: ", e$message)
    })
  }
  
  # ===== PC-Clocks via methylCIPHER =====
  # PC clocks REQUIRE phenotype data with Age and Female columns
  # PC clocks also require external data file from Zenodo
  # Default: Use Horvath2 for Age and InferredSex for Female
  if (requireNamespace("methylCIPHER", quietly = TRUE)) {
    tryCatch({
      if (verbose) message("  Computing PC clocks...")
      
      # Check for and download PC clocks data file if needed
      pc_data_path <- NULL
      cache_dir <- get_manifest_cache_dir()  # Reuse our cache directory
      pc_data_file <- file.path(cache_dir, "PCClocks_data.qs2")
      
      if (file.exists(pc_data_file)) {
        # Check file size - should be ~2GB, reject if too small (partial download)
        file_size <- file.info(pc_data_file)$size
        if (file_size > 1e9) {  # At least 1GB
          pc_data_path <- pc_data_file
          if (verbose) message("    Using cached PC clocks data (", round(file_size/1e9, 2), " GB)")
        } else {
          if (verbose) message("    Cached PC clocks data appears incomplete (", round(file_size/1e6, 1), " MB), re-downloading...")
          unlink(pc_data_file)  # Remove partial file
        }
      }
      
      if (is.null(pc_data_path)) {
        # Try to download from Zenodo - file is ~2GB, needs extended timeout
        if (verbose) message("    Downloading PC clocks data from Zenodo (~2GB, this may take several minutes)...")
        pc_url <- "https://zenodo.org/records/17162604/files/PCClocks_data.qs2?download=1"
        
        # Increase timeout for large file (30 minutes)
        old_timeout <- getOption("timeout")
        options(timeout = 1800)
        
        download_success <- tryCatch({
          download.file(pc_url, pc_data_file, mode = "wb", quiet = !verbose)
          TRUE
        }, error = function(e) {
          if (verbose) message("    Failed to download PC clocks data: ", e$message)
          FALSE
        }, warning = function(w) {
          if (verbose) message("    Download warning: ", w$message)
          TRUE  # Continue anyway, check file after
        })
        
        # Restore original timeout
        options(timeout = old_timeout)
        
        # Verify download succeeded and file is complete
        if (file.exists(pc_data_file)) {
          file_size <- file.info(pc_data_file)$size
          if (file_size > 1e9) {  # At least 1GB
            pc_data_path <- pc_data_file
            if (verbose) message("    PC clocks data downloaded and cached (", round(file_size/1e9, 2), " GB)")
          } else {
            if (verbose) message("    Download incomplete (", round(file_size/1e6, 1), " MB), removing partial file")
            unlink(pc_data_file)
          }
        }
      }
      
      # Load PC clocks CpG data to global env
      tryCatch({
        data("PCClocks_CpGs", package = "methylCIPHER", envir = .GlobalEnv)
      }, error = function(e) NULL, warning = function(w) NULL)
      
      betas_t <- t(betas)  # methylCIPHER wants samples as rows
      sample_ids <- rownames(betas_t)
      
      if (exists("calcPCClocks", envir = asNamespace("methylCIPHER"))) {
        pc_func <- get("calcPCClocks", envir = asNamespace("methylCIPHER"))
        
        # Build pheno data - get Age values
        if (!is.null(pheno) && is.data.frame(pheno) && "Age" %in% colnames(pheno)) {
          age_values <- pheno$Age
          if (verbose) message("    Using provided Age")
        } else if ("Horvath2" %in% colnames(results)) {
          age_values <- results$Horvath2
          if (verbose) message("    Using Horvath2 as Age proxy")
        } else if ("Horvath1" %in% colnames(results)) {
          age_values <- results$Horvath1
          if (verbose) message("    Using Horvath1 as Age proxy")
        } else {
          age_values <- rep(50, length(sample_ids))
          if (verbose) message("    Warning: No age available, using placeholder 50")
        }
        
        # Get Female values
        if (!is.null(pheno) && is.data.frame(pheno) && "Female" %in% colnames(pheno)) {
          female_values <- pheno$Female
          if (verbose) message("    Using provided Female")
        } else if ("InferredSex_Numeric" %in% colnames(results)) {
          female_values <- results$InferredSex_Numeric
          if (verbose) message("    Using InferredSex_Numeric as Female")
        } else {
          female_values <- rep(0.5, length(sample_ids))
          if (verbose) message("    Warning: No sex available, using 0.5 (unknown)")
        }
        
        # Create pheno dataframe - methylCIPHER requires Sample_ID column
        pheno_df <- data.frame(
          Sample_ID = sample_ids,
          Age = age_values,
          Female = female_values,
          stringsAsFactors = FALSE
        )
        rownames(pheno_df) <- sample_ids
        
        if (verbose) {
          message("    PC clocks pheno: Age range = ", 
                  round(min(pheno_df$Age, na.rm=TRUE), 1), " - ", 
                  round(max(pheno_df$Age, na.rm=TRUE), 1),
                  ", Female: ", sum(pheno_df$Female == 1), " F, ",
                  sum(pheno_df$Female == 0), " M")
        }
        
        # Call calcPCClocks with RData path
        pc_result <- NULL
        if (!is.null(pc_data_path)) {
          pc_result <- tryCatch({
            pc_func(betas_t, pheno_df, RData = pc_data_path)
          }, error = function(e) {
            if (verbose) message("    PC clocks failed: ", e$message)
            NULL
          })
        } else {
          if (verbose) message("    PC clocks skipped: data file not available")
        }
        
        if (!is.null(pc_result)) {
          if (verbose) {
            message("    PC result type: ", class(pc_result)[1])
            if (is.data.frame(pc_result)) {
              message("    PC result dims: ", nrow(pc_result), " x ", ncol(pc_result))
              message("    PC result cols: ", paste(head(colnames(pc_result), 10), collapse=", "))
            }
          }
          
          if (is.data.frame(pc_result)) {
            pc_cols <- c("PCHorvath1", "PCHorvath2", "PCHannum", 
                        "PCPhenoAge", "PCGrimAge", "PCDNAmTL")
            added_pc <- c()
            for (col in colnames(pc_result)) {
              if (col %in% pc_cols && is.numeric(pc_result[[col]])) {
                if (nrow(pc_result) == ncol(betas)) {
                  results[[col]] <- pc_result[[col]]
                  added_pc <- c(added_pc, col)
                }
              }
            }
            if (verbose && length(added_pc) > 0) {
              message("    PC clocks added: ", paste(added_pc, collapse = ", "))
            }
          } else if (is.numeric(pc_result)) {
            if (verbose) message("    PC result is numeric vector of length ", length(pc_result))
          }
        } else {
          if (verbose) message("    PC clocks: No result returned")
        }
      }
    }, error = function(e) {
      if (verbose) message("    PC clocks error: ", e$message)
    })
  }
  
  # ===== EpiMitClocks (additional mitotic clocks) =====
  # Only attempt if we don't have epiTOC2 from direct calculation
  epi_pkg <- NULL
  for (pkg_name in c("EpiMitClocks", "epiTOC2")) {
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      epi_pkg <- pkg_name
      break
    }
  }
  
  if (!is.null(epi_pkg) && !"epiTOC2_TNSC" %in% colnames(results)) {
    tryCatch({
      if (verbose) message("  Computing mitotic clocks via EpiMitClocks...")
      
      # Load ALL required data to GLOBAL environment
      epi_data_items <- c("dataETOC3", "cugpmitclockCpG", "epiTOCcpgs3", "estETOC3",
                          "EpiCMITcpgs", "Replitali")
      for (d in epi_data_items) {
        tryCatch({
          data(list = d, package = epi_pkg, envir = .GlobalEnv)
        }, error = function(e) NULL, warning = function(w) NULL)
      }
      
      if (exists("EpiMitClocks", envir = asNamespace(epi_pkg))) {
        tryCatch({
          epi_results <- get("EpiMitClocks", envir = asNamespace(epi_pkg))(data.m = betas, ages.v = NULL)
          
          if (!is.null(epi_results) && is.data.frame(epi_results)) {
            if (nrow(epi_results) == ncol(betas)) {
              for (col in colnames(epi_results)) {
                if (is.numeric(epi_results[[col]]) && !col %in% colnames(results)) {
                  results[[col]] <- epi_results[[col]]
                }
              }
              if (verbose) {
                message("    EpiMitClocks: ", paste(colnames(epi_results), collapse = ", "))
              }
            }
          }
        }, error = function(e) {
          if (verbose) message("    EpiMitClocks failed: ", e$message)
        })
      }
    }, error = function(e) {
      if (verbose) message("    Mitotic clocks error: ", e$message)
    })
  }
  
  # ===== methylCIPHER additional clocks =====
  # Only compute clocks we don't already have from direct calculations
  if (requireNamespace("methylCIPHER", quietly = TRUE)) {
    tryCatch({
      if (verbose) message("  Computing additional clocks via methylCIPHER functions...")
      
      betas_t <- t(betas)  # methylCIPHER wants samples as rows
      
      # Only try clocks we don't have yet
      additional_clocks <- c("AdaptAge", "CausAge", "DamAge", "SystemsAge")
      
      # Load required data
      for (d in c("AdaptAge_CpGs", "CausAge_CpGs", "DamAge_CpGs", "SystemsAge_CpGs")) {
        tryCatch({
          data(list = d, package = "methylCIPHER", envir = .GlobalEnv)
        }, error = function(e) NULL, warning = function(w) NULL)
      }
      
      added_clocks <- c()
      
      for (clock in additional_clocks) {
        if (!clock %in% colnames(results)) {
          func_name <- paste0("calc", clock)
          tryCatch({
            if (exists(func_name, envir = asNamespace("methylCIPHER"))) {
              func <- get(func_name, envir = asNamespace("methylCIPHER"))
              result <- tryCatch({
                func(betas_t, imputation = FALSE)
              }, error = function(e) {
                tryCatch(func(betas_t), error = function(e2) NULL)
              })
              
              if (!is.null(result) && is.numeric(result) && length(result) == ncol(betas)) {
                results[[clock]] <- as.numeric(result)
                added_clocks <- c(added_clocks, clock)
              }
            }
          }, error = function(e) NULL)
        }
      }
      
      if (verbose && length(added_clocks) > 0) {
        message("    Additional clocks: ", paste(added_clocks, collapse = ", "))
      }
    }, error = function(e) {
      if (verbose) message("    Additional clocks error: ", e$message)
    })
  }
  
  return(results)
}