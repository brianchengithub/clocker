# ============================================================================
# Clock Calculation Engines
# Low-level: coefficient handling, weighted-sum calculator, epiTOC2 direct,
# and the missing-probe report writer.
# ============================================================================


#' Horvath age transformation (anti-log transformation)
#'
#' Converts raw Horvath weighted sum to DNAm age in years using the
#' piecewise transformation defined in Horvath (2013).
#'
#' @keywords internal
horvath_age_transform <- function(x, adult_age = 20) {
  ifelse(x < 0,
         (adult_age + 1) * exp(x) - 1,
         (adult_age + 1) * x + adult_age)
}


# ---------------------------------------------------------------------------
# Coverage tracking
#
# Per-clock coverage info is collected into .qc_env$coverage_log so that the
# orchestrator can write a CSV report at the end of the run (feature A2).
# ---------------------------------------------------------------------------

#' Initialize the coverage log for a run
#' @keywords internal
init_coverage_log <- function() {
  .qc_env$coverage_log <- list()
}

#' Record per-sample coverage for one clock
#' @keywords internal
record_coverage <- function(clock_name, sample_ids, n_total,
                              n_present_per_sample) {
  .qc_env$coverage_log[[clock_name]] <- data.frame(
    sample_id        = sample_ids,
    n_cpgs_total     = n_total,
    n_cpgs_present   = n_present_per_sample,
    pct_missing      = round(100 * (n_total - n_present_per_sample) / n_total, 3),
    stringsAsFactors = FALSE
  )
}

#' Get the consolidated coverage log
#' @keywords internal
get_coverage_log <- function() .qc_env$coverage_log %||% list()


# ---------------------------------------------------------------------------
# Weighted-sum clock calculator
# ---------------------------------------------------------------------------

#' Generic weighted-sum clock calculator
#'
#' Computes intercept + sum(beta * weight). Handles the diverse coefficient-
#' table formats used by upstream packages (data frames, matrices, named
#' numeric vectors). Records per-sample CpG coverage in the coverage log.
#'
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coef_df Data frame, matrix, or named vector with CpG/weight info
#' @param transform_func Optional transformation function applied to results
#' @param clock_name Name of clock (used for coverage logging)
#' @return Named vector of clock values, or NULL on failure. Attributes:
#'   `n_cpgs_total`, `n_cpgs_used` (total over all samples).
#' @keywords internal
calc_weighted_sum_clock <- function(betas, coef_df,
                                      transform_func = NULL,
                                      clock_name = NA_character_) {

  if (is.null(coef_df)) return(NULL)

  # Coerce non-data.frame coefficient containers
  if (is.numeric(coef_df) && !is.matrix(coef_df) && !is.data.frame(coef_df)) {
    if (is.null(names(coef_df))) return(NULL)
    coef_df <- data.frame(CpG = names(coef_df),
                          Coefficient = as.numeric(coef_df),
                          stringsAsFactors = FALSE)
  }
  if (is.matrix(coef_df)) {
    coef_df <- as.data.frame(coef_df, stringsAsFactors = FALSE)
    if (!is.null(rownames(coef_df)) && !any(sapply(coef_df, is.character))) {
      coef_df$CpG <- rownames(coef_df)
    }
  }
  if (!is.data.frame(coef_df) || ncol(coef_df) == 0L || nrow(coef_df) == 0L) {
    return(NULL)
  }

  # Auto-detect CpG column
  cpg_candidates <- c("CpG", "CpGmarker", "probe", "Probe", "cpg", "ID",
                       "ProbeID", "probe_id", "CpG_ID", "Marker")
  cpg_col <- intersect(cpg_candidates, colnames(coef_df))[1]
  if (is.na(cpg_col)) {
    char_cols <- which(sapply(coef_df, is.character) | sapply(coef_df, is.factor))
    if (length(char_cols) > 0L) {
      cpg_col <- colnames(coef_df)[char_cols[1]]
    } else if (!is.null(rownames(coef_df)) &&
               any(grepl("^cg|^ch", rownames(coef_df), ignore.case = TRUE))) {
      coef_df$CpG <- rownames(coef_df)
      cpg_col <- "CpG"
    } else {
      cpg_col <- colnames(coef_df)[1]
    }
  }

  # Auto-detect weight column
  weight_candidates <- c("Coefficient", "CoefficientTraining", "weight", "Weight",
                         "coef", "beta", "Beta", "Effect", "effect")
  weight_col <- intersect(weight_candidates, colnames(coef_df))[1]
  if (is.na(weight_col)) {
    numeric_cols <- which(sapply(coef_df, is.numeric))
    numeric_cols <- setdiff(names(numeric_cols), cpg_col)
    if (length(numeric_cols) > 0L) weight_col <- numeric_cols[1]
  }
  if (is.na(weight_col) || is.null(weight_col)) return(NULL)

  cpgs <- as.character(coef_df[[cpg_col]])
  weights <- as.numeric(coef_df[[weight_col]])

  # Intercept handling (FIX M9: warn if multiple intercept rows)
  intercept_mask <- cpgs %in% c("(Intercept)", "Intercept") |
                    is.na(cpgs) | cpgs == ""
  if (any(intercept_mask)) {
    if (sum(intercept_mask) > 1L) {
      warning(sprintf("Clock '%s': %d intercept rows found, using first",
                      clock_name %||% "?", sum(intercept_mask)))
    }
    intercept <- weights[intercept_mask][1]
    if (is.na(intercept)) intercept <- 0
    cpgs <- cpgs[!intercept_mask]
    weights <- weights[!intercept_mask]
  } else {
    intercept <- 0
  }
  if (length(cpgs) == 0L) return(NULL)

  # Match against beta matrix
  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  if (sum(valid) == 0L) return(NULL)

  betas_subset  <- betas[matched_idx[valid], , drop = FALSE]
  weights_valid <- weights[valid]

  # ---- Coverage tracking (per-sample) ----
  n_total <- length(cpgs)
  if (anyNA(betas_subset)) {
    # n_present per sample = number of clock CpGs that are non-NA in that sample
    n_present_per_sample <- colSums(!is.na(betas_subset))
  } else {
    n_present_per_sample <- rep(sum(valid), ncol(betas))
  }
  if (!is.na(clock_name)) {
    record_coverage(clock_name, colnames(betas), n_total, n_present_per_sample)
  }

  # ---- Compute clock value ----
  # Notes on NA handling: imputation runs upstream so betas should be NA-free,
  # but we still apply na.rm=TRUE defensively. If a sample has unexpectedly
  # many NA betas, the coverage log will surface it.
  clock_values <- intercept + colSums(betas_subset * weights_valid, na.rm = TRUE)

  if (!is.null(transform_func)) clock_values <- transform_func(clock_values)

  attr(clock_values, "n_cpgs_total")    <- n_total
  attr(clock_values, "n_cpgs_used_max") <- sum(valid)
  clock_values
}


# ---------------------------------------------------------------------------
# epiTOC2 (mitotic clock)
# ---------------------------------------------------------------------------

#' Calculate epiTOC2 mitotic clock directly
#'
#' Computes total stem cell divisions (TNSC) per the epiTOC2 algorithm.
#' (Per request, this remains as apply()-over-columns rather than full
#' vectorization to keep the inner-loop semantics traceable to the original
#' algorithm description.)
#'
#' @keywords internal
calc_epitoc2_direct <- function(betas, coeffs, verbose = TRUE) {

  estETOC2 <- NULL
  if ("dataETOC3.l" %in% names(coeffs)) {
    data_list <- coeffs[["dataETOC3.l"]]
    if (is.list(data_list) && length(data_list) >= 1L) {
      estETOC2 <- data_list[[1]]
    }
  }
  if (is.null(estETOC2) && "EpiToc2_CpGs" %in% names(coeffs)) {
    estETOC2 <- coeffs$EpiToc2_CpGs
  }
  if (is.null(estETOC2)) return(NULL)

  cpgs <- rownames(estETOC2)
  if (is.null(cpgs)) return(NULL)

  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)
  n_valid <- sum(valid)
  if (n_valid == 0L) return(NULL)

  if (verbose) {
    message(sprintf("    epiTOC2: Using %d of %d CpGs", n_valid, length(cpgs)))
  }

  betas_matched <- betas[matched_idx[valid], , drop = FALSE]
  params_matched <- estETOC2[valid, , drop = FALSE]

  # Coverage tracking
  if (anyNA(betas_matched)) {
    n_present_per_sample <- colSums(!is.na(betas_matched))
  } else {
    n_present_per_sample <- rep(n_valid, ncol(betas))
  }
  record_coverage("epiTOC2_TNSC", colnames(betas), length(cpgs),
                   n_present_per_sample)

  if (ncol(params_matched) >= 2L) {
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
  NULL
}


# ---------------------------------------------------------------------------
# Missing-probe CSV report (feature A1)
# ---------------------------------------------------------------------------

#' Write the per-sample, per-clock missing-probe CSV report
#'
#' One row per (sample, clock); columns include n_total, n_present, pct_missing,
#' plus a wide-format version on the second sheet (samples x clocks of pct_missing).
#' The wide-format version is what most users will inspect.
#'
#' @param path Output CSV path
#' @return Invisibly: the path written
#' @keywords internal
write_missing_probe_report <- function(path, verbose = TRUE) {
  log <- get_coverage_log()
  if (length(log) == 0L) {
    if (verbose) message("    No coverage data available; skipping missing-probe report")
    return(invisible(NULL))
  }

  # Long form
  rows <- do.call(rbind, lapply(names(log), function(clock) {
    df <- log[[clock]]
    df$clock <- clock
    df
  }))
  rows <- rows[, c("sample_id", "clock", "n_cpgs_total",
                   "n_cpgs_present", "pct_missing")]

  # Wide form: samples x clocks, values = pct_missing
  wide <- tryCatch({
    sample_ids <- unique(rows$sample_id)
    clocks <- unique(rows$clock)
    out <- matrix(NA_real_, nrow = length(sample_ids), ncol = length(clocks),
                   dimnames = list(sample_ids, clocks))
    for (cl in clocks) {
      sub <- rows[rows$clock == cl, ]
      out[sub$sample_id, cl] <- sub$pct_missing
    }
    data.frame(sample_id = rownames(out), out, check.names = FALSE,
                stringsAsFactors = FALSE)
  }, error = function(e) NULL)

  utils::write.csv(rows, path, row.names = FALSE)
  if (!is.null(wide)) {
    wide_path <- sub("\\.csv$", "_wide.csv", path, ignore.case = TRUE)
    if (wide_path == path) wide_path <- paste0(path, ".wide.csv")
    utils::write.csv(wide, wide_path, row.names = FALSE)
    if (verbose) {
      log_msg("    Missing-probe report (long):  %s", path,  verbose = TRUE)
      log_msg("    Missing-probe report (wide):  %s", wide_path, verbose = TRUE)
    }
  } else if (verbose) {
    log_msg("    Missing-probe report: %s", path, verbose = TRUE)
  }
  invisible(path)
}
