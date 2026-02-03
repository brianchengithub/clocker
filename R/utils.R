# ============================================================================
# Utility Functions
# Logging, validation, platform detection, core optimization
# ============================================================================


#' Log message with formatting
#' @param fmt Format string (passed to sprintf)
#' @param ... Arguments for format string
#' @param verbose Logical. Print message if TRUE.
#' @keywords internal
log_msg <- function(fmt, ..., verbose = TRUE) {
  if (verbose) {
    message(sprintf(fmt, ...))
  }
}


#' Null coalescing operator
#' @param a First value
#' @param b Fallback value if a is NULL
#' @return a if not NULL, otherwise b
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Validate beta matrix format and values
#' @param betas Matrix to validate
#' @return Validated (possibly corrected) beta matrix
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
    betas <- 2^betas / (2^betas + 1)
  }

  return(betas)
}


#' Detect array platform from probe IDs
#' @param probe_ids Character vector of probe IDs
#' @return Character string identifying the platform
#' @keywords internal
detect_array_platform <- function(probe_ids) {

  n_probes <- length(probe_ids)

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
#' @param n_samples Number of samples
#' @param n_probes Number of probes
#' @param requested_cores User-requested core count (or NULL)
#' @param verbose Print progress
#' @return Integer number of cores to use
#' @keywords internal
determine_optimal_cores <- function(n_samples, n_probes, requested_cores, verbose) {

  available_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(available_cores)) available_cores <- 1

  bytes_per_sample <- n_probes * 8 * 3

  available_ram_gb <- tryCatch({
    if (.Platform$OS.type == "unix") {
      if (file.exists("/proc/meminfo")) {
        meminfo <- readLines("/proc/meminfo", n = 3)
        mem_free <- as.numeric(gsub("[^0-9]", "",
          meminfo[grep("MemAvailable|MemFree", meminfo)[1]]))
        mem_free / 1024 / 1024
      } else {
        mem_str <- system("sysctl -n hw.memsize", intern = TRUE)
        as.numeric(mem_str) / 1024^3
      }
    } else {
      mem_str <- system("wmic OS get FreePhysicalMemory", intern = TRUE)[2]
      as.numeric(trimws(mem_str)) / 1024 / 1024
    }
  }, error = function(e) 8)

  if (is.na(available_ram_gb) || available_ram_gb <= 0) {
    available_ram_gb <- 8
  }

  usable_ram_gb <- max(available_ram_gb - 2, 1)
  gb_per_sample <- bytes_per_sample / 1024^3
  max_parallel_samples <- floor(usable_ram_gb / gb_per_sample)

  memory_limited_cores <- max(1, min(available_cores, floor(max_parallel_samples)))
  sample_limited_cores <- min(available_cores, n_samples)

  optimal_cores <- min(memory_limited_cores, sample_limited_cores, available_cores - 1)
  optimal_cores <- max(1, optimal_cores)

  if (!is.null(requested_cores)) {
    if (requested_cores > optimal_cores && verbose) {
      warning(sprintf(
        "Requested %d cores may exceed available resources. Recommended: %d cores",
        requested_cores, optimal_cores
      ))
    }
    return(max(1, min(requested_cores, available_cores)))
  }

  return(optimal_cores)
}
