# ============================================================================
# Utility Functions
# Logging, validation, platform detection, core optimization, progress bars,
# coercion helpers
# ============================================================================


# ---------------------------------------------------------------------------
# Private package environment (replaces .GlobalEnv writes)
# Initialized in zzz.R but referenced here.
# ---------------------------------------------------------------------------
# .qc_env is created in zzz.R via:  .qc_env <- new.env(parent = emptyenv())
# Functions in this package read/write to .qc_env, never to .GlobalEnv.


#' Log message with formatting
#' @param fmt Format string (passed to sprintf)
#' @param ... Arguments for format string
#' @param verbose Logical. Print message if TRUE.
#' @keywords internal
log_msg <- function(fmt, ..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    message(sprintf(fmt, ...))
  }
}


#' Null coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a


# ---------------------------------------------------------------------------
# Progress bar utilities
# Replaces per-sample message lines with a single updating progress bar.
# Uses base txtProgressBar for portability; optionally uses {progress} for
# nicer output (ETA + custom format) when installed.
# ---------------------------------------------------------------------------

#' Create a progress bar
#'
#' Wraps base::txtProgressBar (or progress::progress_bar if available) into a
#' uniform interface used throughout the package. Returns a list with $tick(n)
#' and $done() closures so callers can be agnostic to the backend.
#'
#' @param total Total number of ticks expected
#' @param label Short label shown to the user (e.g., "Imputing")
#' @param verbose If FALSE, returns no-op closures (silent)
#' @return List with $tick(n=1) and $done() functions
#' @keywords internal
make_progress_bar <- function(total, label = "Progress", verbose = TRUE) {
  if (!isTRUE(verbose) || total < 1) {
    return(list(tick = function(n = 1L) invisible(NULL),
                done = function() invisible(NULL)))
  }

  use_progress_pkg <- requireNamespace("progress", quietly = TRUE) &&
                      interactive()

  if (use_progress_pkg) {
    pb <- progress::progress_bar$new(
      format = sprintf("  %s [:bar] :percent (:current/:total) eta: :eta",
                       label),
      total = total, clear = FALSE, width = 70
    )
    list(
      tick = function(n = 1L) for (i in seq_len(n)) pb$tick(),
      done = function() invisible(NULL)
    )
  } else {
    cat(sprintf("  %s: ", label))
    pb <- utils::txtProgressBar(min = 0, max = total, style = 3, width = 50)
    counter <- 0L
    list(
      tick = function(n = 1L) {
        counter <<- counter + n
        utils::setTxtProgressBar(pb, min(counter, total))
      },
      done = function() {
        close(pb)
        cat("\n")
      }
    )
  }
}


# ---------------------------------------------------------------------------
# Cache directory (CRAN-compliant)
# ---------------------------------------------------------------------------

#' Get the package cache directory (CRAN-compliant)
#'
#' Resolution order (first writable wins):
#'   1. Option `clocker.cache_dir` if set
#'   2. Environment variable `CLOCKER_CACHE` if set
#'   3. tools::R_user_dir("clocker", "cache")  -- respects XDG on Linux
#'   4. tempdir()/clocker  -- always writable, but session-scoped
#'
#' @param subdir Optional subdirectory (e.g., "manifests", "pcclocks")
#' @return Absolute path to writable cache directory
#' @keywords internal
get_cache_dir <- function(subdir = NULL) {

  candidates <- c(
    getOption("clocker.cache_dir"),
    Sys.getenv("CLOCKER_CACHE", unset = NA_character_),
    tryCatch(tools::R_user_dir("clocker", "cache"),
             error = function(e) NA_character_),
    file.path(tempdir(), "clocker")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  for (cand in candidates) {
    target <- if (!is.null(subdir)) file.path(cand, subdir) else cand
    if (dir.exists(target)) return(target)
    ok <- suppressWarnings(
      dir.create(target, recursive = TRUE, showWarnings = FALSE))
    if (ok || dir.exists(target)) return(target)
  }

  # Last resort - tempdir is always writable
  fallback <- file.path(tempdir(), "clocker", subdir %||% "")
  dir.create(fallback, recursive = TRUE, showWarnings = FALSE)
  return(fallback)
}


#' Get user data directory for reference databases
#'
#' Reference databases are user-built and large (>1GB), so they live outside
#' the cache (so cache cleaners don't nuke them).
#'
#' @keywords internal
get_data_dir <- function(subdir = NULL) {
  candidates <- c(
    getOption("clocker.data_dir"),
    Sys.getenv("CLOCKER_DATA", unset = NA_character_),
    tryCatch(tools::R_user_dir("clocker", "data"),
             error = function(e) NA_character_),
    file.path(tempdir(), "clocker_data")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]

  for (cand in candidates) {
    target <- if (!is.null(subdir)) file.path(cand, subdir) else cand
    if (dir.exists(target)) return(target)
    ok <- suppressWarnings(
      dir.create(target, recursive = TRUE, showWarnings = FALSE))
    if (ok || dir.exists(target)) return(target)
  }

  fallback <- file.path(tempdir(), "clocker_data", subdir %||% "")
  dir.create(fallback, recursive = TRUE, showWarnings = FALSE)
  fallback
}


# ---------------------------------------------------------------------------
# Coercion helpers
# ---------------------------------------------------------------------------

#' Coerce a Female/sex column to integer 0/1 (M=0, F=1)
#'
#' Accepts: 0, 1, "M", "F", "Male", "Female" (any case), or any combination.
#' Values like 0.5 ("unknown") are preserved as NA with a warning.
#'
#' @param x Vector of sex/female indicators
#' @param na_on_unknown If TRUE (default), 0.5 / "U" / NA become NA.
#'   If FALSE, they are coerced to 0 (legacy behavior)
#' @return Integer vector of 0L (Male) / 1L (Female) / NA
#' @keywords internal
coerce_female_indicator <- function(x, na_on_unknown = TRUE) {

  if (length(x) == 0) return(integer(0))

  result <- rep(NA_integer_, length(x))

  if (is.numeric(x)) {
    is_zero    <- !is.na(x) & x == 0
    is_one     <- !is.na(x) & x == 1
    is_unknown <- !is.na(x) & !(is_zero | is_one)
    result[is_zero] <- 0L
    result[is_one]  <- 1L
    if (any(is_unknown) && !na_on_unknown) result[is_unknown] <- 0L
    if (any(is_unknown) && na_on_unknown) {
      warning(sprintf(
        "%d Female value(s) outside {0, 1} (e.g., 0.5) coerced to NA",
        sum(is_unknown)))
    }
  } else {
    x_chr <- tolower(trimws(as.character(x)))
    male_set    <- c("0", "m", "male")
    female_set  <- c("1", "f", "female")
    unknown_set <- c("", "u", "unknown", "na", "0.5", "5e-1")

    result[x_chr %in% male_set]   <- 0L
    result[x_chr %in% female_set] <- 1L

    bad <- !is.na(x_chr) & !(x_chr %in% c(male_set, female_set, unknown_set))
    if (any(bad)) {
      warning(sprintf(
        "Unrecognized sex value(s) coerced to NA: %s",
        paste(unique(x[bad]), collapse = ", ")))
    }
  }

  return(result)
}


# ---------------------------------------------------------------------------
# Beta matrix validation
# ---------------------------------------------------------------------------

#' Validate beta matrix format and values
#'
#' Robust M-value detection: only flips M->beta when the *median* (not the
#' min/max) of the input falls outside the unit interval (0 to 1). A
#' single bad probe in 800k shouldn't wreck an entire matrix.
#'
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

  # Probe-name sniff
  probe_pattern <- sum(grepl("^cg|^ch", rownames(betas), ignore.case = TRUE))
  if (probe_pattern < nrow(betas) * 0.5) {
    warning("Less than 50% of rownames appear to be CpG probe IDs. ",
            "Expected format: 'cg00000029', 'ch.1.1234', etc.")
  }

  # ---- Robust M-value vs beta detection ----
  # Use median (and IQR midpoint) instead of min/max; extreme outliers in
  # 800k probes can spuriously trigger the old min/max test.
  finite_vals <- betas[is.finite(betas)]
  if (length(finite_vals) > 0) {
    val_median <- stats::median(finite_vals)
    val_q1 <- stats::quantile(finite_vals, 0.25, names = FALSE)
    val_q3 <- stats::quantile(finite_vals, 0.75, names = FALSE)

    looks_like_mvals <- (val_median < -0.5 || val_median > 1.5) ||
                        (val_q1 < -0.2 && val_q3 > 1.2)

    if (looks_like_mvals) {
      message("Input median = ", round(val_median, 3),
              " (IQR ", round(val_q1, 2), "-", round(val_q3, 2), ").")
      message("Values appear to be M-values; converting via 2^x / (2^x + 1).")
      betas <- 2^betas / (2^betas + 1)
    } else if (min(finite_vals) < -0.05 || max(finite_vals) > 1.05) {
      # Beta-like overall but with stray out-of-range values: clip
      n_low  <- sum(finite_vals < 0)
      n_high <- sum(finite_vals > 1)
      if (n_low + n_high > 0) {
        message(sprintf(
          "Clipping %d out-of-range beta values (%d < 0, %d > 1) to [0, 1]",
          n_low + n_high, n_low, n_high))
        betas[!is.na(betas) & betas < 0] <- 0
        betas[!is.na(betas) & betas > 1] <- 1
      }
    }
  }

  return(betas)
}


# ---------------------------------------------------------------------------
# EPICv2 normalization
# ---------------------------------------------------------------------------

#' Normalize EPICv2 probe names to base CpG IDs
#'
#' EPICv2 arrays use probe names with strand/channel/replicate suffixes
#' (e.g., cg00000029_TC11). Most clock algorithms, EpiDISH references, and
#' DunedinPACE expect base CpG IDs (e.g., cg00000029).
#'
#' For CpGs with multiple replicates, this function tries to select the
#' replicate whose genomic address matches the original EPIC v1 or 450K
#' probe -- preserving cross-platform comparability. Replicates that cannot
#' be disambiguated by position are averaged.
#'
#' @param betas Beta matrix (probes as rows, samples as columns)
#' @param verbose Print progress
#' @return Beta matrix with base CpG IDs as rownames
#' @keywords internal
normalize_epicv2_probes <- function(betas, verbose = TRUE) {

  probe_names <- rownames(betas)
  has_suffix <- grepl("_[TB]C[0-9]{2}$", probe_names)

  # FIX C3: don't gate on a 10% threshold. Even a single suffixed probe
  # means we must normalize, otherwise that probe's value is silently
  # dropped from clock matching. Performance cost on tiny suffix sets is
  # negligible because the function returns early below if there are no
  # multi-probe collisions.
  if (!any(has_suffix)) {
    return(betas)
  }

  if (verbose) message("  Normalizing EPICv2 probe names to base CpG IDs...")

  base_names <- sub("_[TB]C[0-9]{2}$", "", probe_names)

  base_counts  <- table(base_names)
  multi_bases  <- names(base_counts[base_counts > 1])
  single_bases <- names(base_counts[base_counts == 1])
  unique_bases <- names(base_counts)

  if (verbose) {
    message(sprintf("    %d probes -> %d unique base CpGs",
                    length(probe_names), length(unique_bases)))
    message(sprintf("    %d CpGs with replicate probes", length(multi_bases)))
  }

  if (length(multi_bases) == 0) {
    rownames(betas) <- base_names
    return(betas)
  }

  position_resolved <- resolve_replicates_by_position(
    probe_names, base_names, multi_bases, verbose
  )

  new_betas <- matrix(NA_real_,
                      nrow = length(unique_bases),
                      ncol = ncol(betas),
                      dimnames = list(unique_bases, colnames(betas)))

  # Single-probe bases: direct copy
  single_idx <- which(base_names %in% single_bases)
  new_betas[base_names[single_idx], ] <- betas[single_idx, , drop = FALSE]

  n_pos_resolved <- 0L
  n_averaged <- 0L

  for (base in multi_bases) {
    idx <- which(base_names == base)

    if (base %in% names(position_resolved)) {
      best_probe <- position_resolved[[base]]
      best_idx <- which(probe_names == best_probe)
      if (length(best_idx) == 1) {
        new_betas[base, ] <- betas[best_idx, , drop = FALSE]
        n_pos_resolved <- n_pos_resolved + 1L
        next
      }
    }

    # Fallback: average replicates
    if (ncol(betas) == 1) {
      new_betas[base, 1] <- mean(betas[idx, 1], na.rm = TRUE)
    } else {
      new_betas[base, ] <- colMeans(betas[idx, , drop = FALSE], na.rm = TRUE)
    }
    n_averaged <- n_averaged + 1L
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
#' Matching priority:
#'   1. Probe_beg + DESIGN + probe_strand
#'   2. Probe_beg + DESIGN
#'   3. Probe_beg only
#'   4. DESIGN + probe_strand
#'   5. No match -> caller averages
#'
#' @keywords internal
resolve_replicates_by_position <- function(probe_names, base_names,
                                            multi_bases, verbose) {

  resolved <- list()

  epicv2_manifest <- tryCatch(download_manifest("EPICv2", verbose = FALSE),
                               error = function(e) NULL)
  if (is.null(epicv2_manifest)) {
    if (verbose) message("    Could not load EPICv2 manifest for matching")
    return(resolved)
  }

  legacy_epic <- tryCatch(download_manifest("EPIC", verbose = FALSE),
                           error = function(e) NULL)
  legacy_450k <- tryCatch(download_manifest("HM450", verbose = FALSE),
                           error = function(e) NULL)

  if (is.null(legacy_epic) && is.null(legacy_450k)) {
    if (verbose) message("    Could not load EPIC or 450K manifests for matching")
    return(resolved)
  }

  find_col <- function(df, candidates) {
    for (col in candidates) if (col %in% colnames(df)) return(col)
    return(NULL)
  }

  probe_id_names  <- c("Probe_ID", "probeID", "IlmnID", "Name")
  probe_beg_names <- c("Probe_beg", "probe_beg", "Probe_Start", "PROBE_BEG")
  design_names    <- c("DESIGN", "design", "Infinium_Design_Type", "Type")
  strand_names    <- c("probe_strand", "Strand", "strand", "STRAND")

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

    info <- data.frame(probe = target_probes[valid],
                       stringsAsFactors = FALSE)
    info$probe_beg <- if (!is.null(b_col)) manifest[[b_col]][idx[valid]] else NA
    info$design    <- if (!is.null(d_col)) manifest[[d_col]][idx[valid]] else NA
    info$strand    <- if (!is.null(s_col)) manifest[[s_col]][idx[valid]] else NA
    info
  }

  multi_full_probes <- probe_names[base_names %in% multi_bases]
  v2_info <- extract_design_info(epicv2_manifest, multi_full_probes)
  if (is.null(v2_info) || nrow(v2_info) == 0) {
    if (verbose) message("    No EPICv2 replicate probes found in manifest")
    return(resolved)
  }
  v2_info$base <- sub("_[TB]C[0-9]{2}$", "", v2_info$probe)

  legacy_info <- NULL
  if (!is.null(legacy_epic)) {
    legacy_info <- extract_design_info(legacy_epic, multi_bases)
  }
  if (!is.null(legacy_450k)) {
    covered <- if (!is.null(legacy_info)) legacy_info$probe else character(0)
    uncovered <- setdiff(multi_bases, covered)
    if (length(uncovered) > 0) {
      info_450k <- extract_design_info(legacy_450k, uncovered)
      if (!is.null(info_450k)) {
        legacy_info <- if (is.null(legacy_info)) info_450k
                       else rbind(legacy_info, info_450k)
      }
    }
  }
  if (is.null(legacy_info) || nrow(legacy_info) == 0) return(resolved)

  rownames(legacy_info) <- legacy_info$probe

  for (base in multi_bases) {
    if (!base %in% legacy_info$probe) next
    leg <- legacy_info[base, ]
    reps <- v2_info[v2_info$base == base, , drop = FALSE]
    if (nrow(reps) == 0) next

    best <- NULL
    if (!is.na(leg$probe_beg) && !is.na(leg$design) && !is.na(leg$strand)) {
      hits <- which(reps$probe_beg == leg$probe_beg &
                    reps$design == leg$design &
                    reps$strand == leg$strand)
      if (length(hits) == 1) best <- reps$probe[hits[1]]
    }
    if (is.null(best) && !is.na(leg$probe_beg) && !is.na(leg$design)) {
      hits <- which(reps$probe_beg == leg$probe_beg &
                    reps$design == leg$design)
      if (length(hits) >= 1) best <- reps$probe[hits[1]]
    }
    if (is.null(best) && !is.na(leg$probe_beg)) {
      hits <- which(reps$probe_beg == leg$probe_beg)
      if (length(hits) >= 1) best <- reps$probe[hits[1]]
    }
    if (is.null(best) && !is.na(leg$design) && !is.na(leg$strand)) {
      hits <- which(reps$design == leg$design & reps$strand == leg$strand)
      if (length(hits) == 1) best <- reps$probe[hits[1]]
    }

    if (!is.null(best)) resolved[[base]] <- best
  }

  resolved
}


# ---------------------------------------------------------------------------
# Platform detection
# ---------------------------------------------------------------------------

#' Detect array platform from probe IDs
#' @keywords internal
detect_array_platform <- function(probe_ids) {

  n_probes <- length(probe_ids)

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


# ---------------------------------------------------------------------------
# Resource estimation
# ---------------------------------------------------------------------------

#' Determine optimal number of cores
#'
#' The "* 3" multiplier on per-sample memory accounts for: (1) input matrix,
#' (2) intermediate copy used during clock projection, (3) scratch space for
#' matrix-multiply temporaries.
#'
#' @keywords internal
determine_optimal_cores <- function(n_samples, n_probes,
                                     requested_cores, verbose) {

  available_cores <- parallel::detectCores(logical = FALSE)
  if (is.na(available_cores)) available_cores <- 1L

  # Memory estimate: n_probes * 8 bytes/double * 3 (input + intermediate + scratch)
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

  if (is.na(available_ram_gb) || available_ram_gb <= 0) available_ram_gb <- 8

  usable_ram_gb <- max(available_ram_gb - 2, 1)
  gb_per_sample <- bytes_per_sample / 1024^3
  max_parallel <- floor(usable_ram_gb / gb_per_sample)

  memory_limited_cores <- max(1L, min(available_cores, floor(max_parallel)))
  sample_limited_cores <- min(available_cores, n_samples)

  optimal <- min(memory_limited_cores, sample_limited_cores, available_cores - 1L)
  optimal <- max(1L, as.integer(optimal))

  if (!is.null(requested_cores)) {
    if (requested_cores > optimal && verbose) {
      warning(sprintf(
        "Requested %d cores may exceed available resources. Recommended: %d",
        requested_cores, optimal))
    }
    return(max(1L, min(as.integer(requested_cores), available_cores)))
  }

  optimal
}


# ---------------------------------------------------------------------------
# File-format helpers (qs2-aware I/O)
# ---------------------------------------------------------------------------

#' Save an R object using qs2 if available, fall back to RDS
#' @keywords internal
qc_save <- function(object, path) {
  if (requireNamespace("qs2", quietly = TRUE) && grepl("\\.qs2?$", path)) {
    qs2::qs_save(object, path)
  } else {
    if (grepl("\\.qs2?$", path)) path <- sub("\\.qs2?$", ".rds", path)
    saveRDS(object, path)
  }
  invisible(path)
}

#' Load an R object using qs2 if available, fall back to RDS
#' @keywords internal
qc_load <- function(path) {
  if (grepl("\\.qs2?$", path) && requireNamespace("qs2", quietly = TRUE)) {
    return(qs2::qs_read(path))
  }
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  if (file.exists(sub("\\.qs2?$", ".rds", path))) {
    return(readRDS(sub("\\.qs2?$", ".rds", path)))
  }
  stop("Cannot load: ", path)
}


#' SHA-256 hash of a file (returns NA if 'digest' not installed)
#' @keywords internal
file_sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}


#' Download a URL to a file with retries and an optional SHA-256 check
#' @keywords internal
download_with_retry <- function(url, dest, retries = 2L,
                                 expected_sha256 = NULL,
                                 verbose = TRUE) {
  use_curl <- requireNamespace("curl", quietly = TRUE)

  for (attempt in 0:retries) {
    res <- tryCatch({
      if (use_curl) {
        curl::curl_download(url, dest, mode = "wb", quiet = !verbose)
      } else {
        utils::download.file(url, dest, mode = "wb", quiet = !verbose)
      }
      TRUE
    }, error = function(e) e, warning = function(w) w)

    if (isTRUE(res) && file.exists(dest) && file.info(dest)$size > 0) {
      if (!is.null(expected_sha256) && !is.na(expected_sha256)) {
        got <- file_sha256(dest)
        if (!is.na(got) && !identical(got, expected_sha256)) {
          if (verbose) {
            message("    SHA-256 mismatch: expected ", expected_sha256,
                    ", got ", got, " - re-trying")
          }
          unlink(dest)
        } else {
          return(invisible(TRUE))
        }
      } else {
        return(invisible(TRUE))
      }
    }
    if (attempt < retries) Sys.sleep(2 ^ attempt)
  }
  stop("Failed to download after ", retries + 1L, " attempts: ", url)
}
