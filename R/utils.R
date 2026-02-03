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


#' Normalize EPICv2 probe names to base CpG IDs
#'
#' EPICv2 arrays use probe names with strand/channel/replicate suffixes
#' (e.g., cg00000029_TC11). Most clock algorithms, EpiDISH references, and
#' DunedinPACE expect base CpG IDs (e.g., cg00000029).
#'
#' For CpGs with multiple replicates on EPICv2, this function first tries to
#' select the replicate whose genomic address (chromosome + position) matches
#' the original EPIC v1 or 450K probe — preserving cross-platform comparability.
#' Only replicates that cannot be disambiguated by position are averaged.
#'
#' @param betas Beta matrix (probes as rows, samples as columns)
#' @param verbose Print progress
#' @return Beta matrix with base CpG IDs as rownames (fewer rows if replicates merged)
#' @keywords internal
normalize_epicv2_probes <- function(betas, verbose = TRUE) {

  probe_names <- rownames(betas)

  # EPICv2 suffixes: _TC11, _TC12, _TC21, _TC22, _BC11, _BC12, _BC21, _BC22
  has_suffix <- grepl("_[TB]C[0-9]{2}$", probe_names)

  if (sum(has_suffix) < length(probe_names) * 0.1) {
    return(betas)
  }

  if (verbose) message("  Normalizing EPICv2 probe names to base CpG IDs...")

  base_names <- sub("_[TB]C[0-9]{2}$", "", probe_names)

  base_counts <- table(base_names)
  multi_bases <- names(base_counts[base_counts > 1])
  single_bases <- names(base_counts[base_counts == 1])
  unique_bases <- names(base_counts)

  if (verbose) {
    message(sprintf("    %d EPICv2 probes -> %d unique base CpGs",
                    length(probe_names), length(unique_bases)))
    message(sprintf("    %d CpGs with replicate probes", length(multi_bases)))
  }

  if (length(multi_bases) == 0) {
    rownames(betas) <- base_names
    return(betas)
  }

  # ==========================================================================
  # Position-based replicate selection for multi-probe CpGs
  # Try to match each EPICv2 replicate to the original EPIC/450K address
  # ==========================================================================

  position_resolved <- resolve_replicates_by_position(
    probe_names, base_names, multi_bases, verbose
  )

  # Build new matrix
  new_betas <- matrix(NA_real_, nrow = length(unique_bases), ncol = ncol(betas),
                      dimnames = list(unique_bases, colnames(betas)))

  # Single-probe bases: direct copy
  single_idx <- which(base_names %in% single_bases)
  new_betas[base_names[single_idx], ] <- betas[single_idx, , drop = FALSE]

  # Multi-probe bases
  n_pos_resolved <- 0
  n_averaged <- 0

  for (base in multi_bases) {
    idx <- which(base_names == base)

    if (base %in% names(position_resolved)) {
      # Use the specific replicate that matches the legacy address
      best_probe <- position_resolved[[base]]
      best_idx <- which(probe_names == best_probe)
      if (length(best_idx) == 1) {
        new_betas[base, ] <- betas[best_idx, , drop = FALSE]
        n_pos_resolved <- n_pos_resolved + 1
        next
      }
    }

    # Fallback: average replicates
    if (ncol(betas) == 1) {
      new_betas[base, 1] <- mean(betas[idx, 1], na.rm = TRUE)
    } else {
      new_betas[base, ] <- colMeans(betas[idx, , drop = FALSE], na.rm = TRUE)
    }
    n_averaged <- n_averaged + 1
  }

  if (verbose) {
    message(sprintf("    Resolved %d CpGs by position match to EPIC/450K",
                    n_pos_resolved))
    message(sprintf("    Averaged %d CpGs with no position match", n_averaged))
  }

  return(new_betas)
}


#' Resolve EPICv2 replicate probes by matching probe design to EPIC/450K
#'
#' EPICv2 replicate probes for the same CpG may use different probe designs
#' (different hybridization positions, Infinium Type I vs II, or strand).
#' This function identifies which replicate matches the original EPIC v1 or
#' 450K probe design by comparing Probe_beg (probe hybridization start),
#' DESIGN type, and probe_strand from the manifests.
#'
#' Matching priority:
#'   1. Probe_beg + DESIGN + probe_strand (exact probe design match)
#'   2. Probe_beg + DESIGN (same position and chemistry)
#'   3. Probe_beg only (same hybridization site)
#'   4. DESIGN + probe_strand (same chemistry, different position)
#'   5. No match -> will be averaged by caller
#'
#' @param probe_names Full EPICv2 probe names (with suffixes)
#' @param base_names Base CpG names (suffixes stripped)
#' @param multi_bases Character vector of base CpGs that have multiple replicates
#' @param verbose Print progress
#' @return Named list: base_cpg -> best EPICv2 probe name (full, with suffix)
#' @keywords internal
resolve_replicates_by_position <- function(probe_names, base_names, multi_bases, verbose) {

  resolved <- list()

  # ---- Load manifests ----
  epicv2_manifest <- tryCatch(
    download_manifest("EPICv2", verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(epicv2_manifest)) {
    if (verbose) message("    Could not load EPICv2 manifest for probe matching")
    return(resolved)
  }

  legacy_epic <- tryCatch(download_manifest("EPIC", verbose = FALSE),
                           error = function(e) NULL)
  legacy_450k <- tryCatch(download_manifest("HM450", verbose = FALSE),
                           error = function(e) NULL)

  if (is.null(legacy_epic) && is.null(legacy_450k)) {
    if (verbose) message("    Could not load EPIC or 450K manifests for probe matching")
    return(resolved)
  }

  # ---- Helper: find column by candidate names ----
  find_col <- function(df, candidates) {
    for (col in candidates) {
      if (col %in% colnames(df)) return(col)
    }
    return(NULL)
  }

  probe_id_names <- c("Probe_ID", "probeID", "IlmnID", "Name")
  probe_beg_names <- c("Probe_beg", "probe_beg", "Probe_Start", "PROBE_BEG")
  design_names <- c("DESIGN", "design", "Infinium_Design_Type", "Type")
  strand_names <- c("probe_strand", "Strand", "strand", "STRAND")

  # ---- Extract probe design info from a manifest (vectorized) ----
  extract_design_info <- function(manifest, target_probes) {
    p_col <- find_col(manifest, probe_id_names)
    b_col <- find_col(manifest, probe_beg_names)
    d_col <- find_col(manifest, design_names)
    s_col <- find_col(manifest, strand_names)

    if (is.null(p_col)) return(NULL)

    probes <- manifest[[p_col]]
    idx <- match(target_probes, probes)
    valid <- !is.na(idx)

    if (sum(valid) == 0) return(NULL)

    info <- data.frame(
      probe = target_probes[valid],
      stringsAsFactors = FALSE
    )

    info$probe_beg <- if (!is.null(b_col)) manifest[[b_col]][idx[valid]] else NA
    info$design    <- if (!is.null(d_col)) manifest[[d_col]][idx[valid]] else NA
    info$strand    <- if (!is.null(s_col)) manifest[[s_col]][idx[valid]] else NA

    return(info)
  }

  # ---- EPICv2 replicate probe info ----
  multi_full_probes <- probe_names[base_names %in% multi_bases]
  v2_info <- extract_design_info(epicv2_manifest, multi_full_probes)

  if (is.null(v2_info) || nrow(v2_info) == 0) {
    if (verbose) message("    No EPICv2 replicate probes found in manifest")
    return(resolved)
  }

  # Add base name column
  v2_info$base <- sub("_[TB]C[0-9]{2}$", "", v2_info$probe)

  # ---- Legacy probe info (EPIC preferred, 450K fallback) ----
  legacy_info <- NULL
  if (!is.null(legacy_epic)) {
    legacy_info <- extract_design_info(legacy_epic, multi_bases)
  }

  # Add 450K probes not already covered by EPIC
  if (!is.null(legacy_450k)) {
    covered <- if (!is.null(legacy_info)) legacy_info$probe else character(0)
    uncovered <- setdiff(multi_bases, covered)
    if (length(uncovered) > 0) {
      info_450k <- extract_design_info(legacy_450k, uncovered)
      if (!is.null(info_450k)) {
        legacy_info <- if (is.null(legacy_info)) info_450k else rbind(legacy_info, info_450k)
      }
    }
  }

  if (is.null(legacy_info) || nrow(legacy_info) == 0) {
    if (verbose) message("    No legacy probe info found for replicate resolution")
    return(resolved)
  }

  # ---- Match: for each multi-CpG, find the best EPICv2 replicate ----
  # Build a lookup for legacy probes keyed by base CpG name
  rownames(legacy_info) <- legacy_info$probe

  for (base in multi_bases) {
    if (!base %in% legacy_info$probe) next

    leg <- legacy_info[base, ]
    reps <- v2_info[v2_info$base == base, , drop = FALSE]

    if (nrow(reps) == 0) next

    best <- NULL

    # Priority 1: Probe_beg + DESIGN + strand all match
    if (!is.na(leg$probe_beg) && !is.na(leg$design) && !is.na(leg$strand)) {
      hits <- which(reps$probe_beg == leg$probe_beg &
                    reps$design == leg$design &
                    reps$strand == leg$strand)
      if (length(hits) == 1) best <- reps$probe[hits[1]]
    }

    # Priority 2: Probe_beg + DESIGN
    if (is.null(best) && !is.na(leg$probe_beg) && !is.na(leg$design)) {
      hits <- which(reps$probe_beg == leg$probe_beg &
                    reps$design == leg$design)
      if (length(hits) >= 1) best <- reps$probe[hits[1]]
    }

    # Priority 3: Probe_beg only (same hybridization site)
    if (is.null(best) && !is.na(leg$probe_beg)) {
      hits <- which(reps$probe_beg == leg$probe_beg)
      if (length(hits) >= 1) best <- reps$probe[hits[1]]
    }

    # Priority 4: DESIGN + strand (same chemistry, different position)
    if (is.null(best) && !is.na(leg$design) && !is.na(leg$strand)) {
      hits <- which(reps$design == leg$design & reps$strand == leg$strand)
      if (length(hits) == 1) best <- reps$probe[hits[1]]
    }

    if (!is.null(best)) {
      resolved[[base]] <- best
    }
  }

  if (verbose) {
    has_beg <- !is.na(legacy_info$probe_beg[1])
    has_des <- !is.na(legacy_info$design[1])
    has_str <- !is.na(legacy_info$strand[1])
    msg_parts <- c()
    if (has_beg) msg_parts <- c(msg_parts, "Probe_beg")
    if (has_des) msg_parts <- c(msg_parts, "DESIGN")
    if (has_str) msg_parts <- c(msg_parts, "strand")
    if (length(msg_parts) > 0) {
      message(sprintf("    Matching criteria: %s", paste(msg_parts, collapse = " + ")))
    }
  }

  return(resolved)
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
