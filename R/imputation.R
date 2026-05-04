# ============================================================================
# kNN Imputation against reference_betas.rds
#
# Reference format (your existing file):
#   reference_betas.rds  -- a NAMED NUMERIC VECTOR of mean beta values,
#                           keyed by CpG ID (length ~ 850k).
#
# Two modes:
#
#   (1) Standard kNN within the target batch
#       For each sample with missing values, find the k most similar OTHER
#       samples in the same input batch (Pearson correlation over probes
#       both samples have non-missing) and impute each missing probe as a
#       Gaussian-kernel-weighted average from those neighbours' values.
#       Missing probes that the neighbours also lack (or for single-sample
#       runs where there are no neighbours) fall through to (2).
#
#   (2) Zero-shot from reference probe means
#       For samples with > zero_shot_threshold fraction missing, OR for
#       residual missing probes that kNN cannot fill, the per-probe mean
#       from reference_betas.rds is used.
#
# This is intentionally simpler than a full reference-database workflow.
# When you later build a multi-sample reference with metadata, the
# imputation engine can be swapped to a "kNN against reference samples"
# variant without touching the rest of the pipeline.
# ============================================================================


REFERENCE_BETAS_FILE <- "reference_betas.rds"


# ---------------------------------------------------------------------------
# Reference-mean lookup
# ---------------------------------------------------------------------------

#' Resolve the reference_betas.rds path
#'
#' Resolution order:
#'   1. Explicit path argument
#'   2. Option `clocker.reference_betas`
#'   3. Env var `CLOCKER_REFERENCE_BETAS`
#'   4. inst/extdata/reference_betas.rds within the installed package
#'   5. data/ directory of the installed package
#' @keywords internal
resolve_reference_betas_path <- function(path = NULL) {
  if (!is.null(path) && nzchar(path)) {
    if (file.exists(path)) return(path)
    stop("reference_betas file not found: ", path)
  }

  cand <- c(
    getOption("clocker.reference_betas"),
    Sys.getenv("CLOCKER_REFERENCE_BETAS", unset = NA_character_),
    tryCatch(system.file("extdata", REFERENCE_BETAS_FILE, package = "clocker"),
              error = function(e) ""),
    tryCatch(system.file("data", REFERENCE_BETAS_FILE, package = "clocker"),
              error = function(e) "")
  )
  cand <- cand[!is.na(cand) & nzchar(cand)]
  for (c in cand) if (file.exists(c)) return(c)
  return(NA_character_)
}


#' Load reference_betas.rds and validate format
#'
#' Expected: a named numeric vector with CpG IDs as names. For backward
#' compatibility we also accept a 1-column data frame, a 1-column matrix,
#' or a list with a `means` element.
#'
#' @keywords internal
load_reference_betas <- function(path = NULL, verbose = TRUE) {

  resolved <- resolve_reference_betas_path(path)
  if (is.na(resolved) || !file.exists(resolved)) {
    if (verbose) {
      message("    reference_betas.rds not found. Set its location via:")
      message("      options(clocker.reference_betas = '/path/to/reference_betas.rds')")
      message("    or place it in inst/extdata/ within the installed package.")
    }
    return(NULL)
  }

  obj <- tryCatch(readRDS(resolved), error = function(e) {
    if (verbose) message("    Failed to read reference_betas: ", e$message)
    NULL
  })
  if (is.null(obj)) return(NULL)

  ref <- NULL
  if (is.numeric(obj) && !is.null(names(obj))) {
    ref <- obj
  } else if (is.data.frame(obj) && nrow(obj) > 0L) {
    num_cols <- which(vapply(obj, is.numeric, logical(1)))
    if (length(num_cols) >= 1L) {
      ref <- obj[[num_cols[1]]]
      names(ref) <- if (!is.null(rownames(obj))) rownames(obj)
                    else as.character(obj[[1]])
    }
  } else if (is.matrix(obj) && ncol(obj) >= 1L && !is.null(rownames(obj))) {
    ref <- obj[, 1]
    names(ref) <- rownames(obj)
  } else if (is.list(obj) && !is.null(obj$means)) {
    ref <- obj$means
  }

  if (is.null(ref) || is.null(names(ref)) || length(ref) == 0L) {
    if (verbose) {
      message("    reference_betas.rds has unrecognized format. ",
              "Expected a named numeric vector of probe means.")
    }
    return(NULL)
  }

  if (verbose) {
    log_msg("    Loaded reference_betas: %d probes from %s",
            length(ref), resolved, verbose = TRUE)
  }
  ref
}


# ---------------------------------------------------------------------------
# Main imputation entry
# ---------------------------------------------------------------------------

#' kNN imputation with zero-shot fallback to reference probe means
#'
#' @param betas Beta matrix (probes x samples). Missing values = NA.
#' @param reference_path Optional path to reference_betas.rds. If NULL,
#'   resolved by `resolve_reference_betas_path()`.
#' @param k Number of nearest in-batch neighbours used for kNN.
#' @param zero_shot_threshold Per-sample fraction of probes missing above
#'   which kNN is bypassed and reference means are used directly.
#' @param distance Distance metric: "pearson" (1 - r) or "euclidean".
#' @param verbose Print progress
#' @return List with `betas` (imputed) and `sample_info` (per-sample frame).
#' @keywords internal
knn_impute <- function(betas,
                        reference_path     = NULL,
                        k                   = 10L,
                        zero_shot_threshold = 0.10,
                        distance            = c("pearson", "euclidean"),
                        verbose             = TRUE) {

  distance <- match.arg(distance)
  k <- as.integer(k)
  if (k < 1L) stop("`k` must be >= 1")
  if (zero_shot_threshold < 0 || zero_shot_threshold > 1) {
    stop("`zero_shot_threshold` must be in [0, 1]")
  }

  n_samples <- ncol(betas)
  n_probes  <- nrow(betas)

  diag_df <- data.frame(
    sample_id            = colnames(betas),
    n_probes_total       = n_probes,
    n_missing_input      = colSums(is.na(betas)),
    pct_missing_input    = round(100 * colSums(is.na(betas)) / n_probes, 3),
    mode                 = NA_character_,
    n_neighbors_used     = NA_integer_,
    mean_neighbor_dist   = NA_real_,
    n_imputed_knn        = 0L,
    n_imputed_zeroshot   = 0L,
    stringsAsFactors     = FALSE
  )

  if (sum(is.na(betas)) == 0L) {
    if (verbose) log_msg("    No missing values to impute", verbose = verbose)
    diag_df$mode <- "none"
    return(list(betas = betas, sample_info = diag_df))
  }

  ref_means <- load_reference_betas(reference_path, verbose = verbose)
  if (is.null(ref_means) && n_samples < 2L) {
    if (verbose) {
      message("    Cannot impute: no reference_betas.rds and only 1 sample. ",
              "Filling missing values with 0.5.")
    }
    betas[is.na(betas)] <- 0.5
    diag_df$mode <- "fallback_constant"
    diag_df$n_imputed_zeroshot <- diag_df$n_missing_input
    return(list(betas = betas, sample_info = diag_df))
  }

  # ---- Per-sample loop with progress bar ----
  pb <- make_progress_bar(n_samples, label = "Imputing", verbose = verbose)
  on.exit(pb$done(), add = TRUE)

  for (i in seq_len(n_samples)) {
    sample_vec <- betas[, i]
    missing_mask <- is.na(sample_vec)
    if (!any(missing_mask)) {
      diag_df$mode[i] <- "none"
      pb$tick()
      next
    }

    pct_missing <- mean(missing_mask)

    # ---- Zero-shot path ----
    if (pct_missing > zero_shot_threshold || n_samples < 2L) {
      diag_df$mode[i] <- "zero_shot"
      filled <- impute_from_reference(sample_vec, missing_mask, ref_means)
      betas[, i] <- filled$values
      diag_df$n_imputed_zeroshot[i] <- filled$n_filled
      pb$tick()
      next
    }

    # ---- kNN within-batch ----
    diag_df$mode[i] <- "knn"
    nn <- find_nearest_neighbors(betas, target_idx = i, k = k,
                                    distance = distance)
    if (is.null(nn) || length(nn$idx) == 0L) {
      # No usable neighbours -- fall back to zero-shot
      diag_df$mode[i] <- "knn_failed_zeroshot"
      filled <- impute_from_reference(sample_vec, missing_mask, ref_means)
      betas[, i] <- filled$values
      diag_df$n_imputed_zeroshot[i] <- filled$n_filled
      pb$tick()
      next
    }

    diag_df$n_neighbors_used[i]   <- length(nn$idx)
    diag_df$mean_neighbor_dist[i] <- mean(nn$dist, na.rm = TRUE)

    # Gaussian-kernel weights: bandwidth = median neighbour distance
    bw <- max(stats::median(nn$dist, na.rm = TRUE), .Machine$double.eps)
    weights <- exp(-(nn$dist / bw)^2)
    weights <- weights / sum(weights)

    # Imputation: weighted mean across neighbours, NA-aware
    nn_block <- betas[missing_mask, nn$idx, drop = FALSE]
    imputed <- weighted_rowmean_with_na(nn_block, weights)

    # Probes still NA after kNN (all neighbours missing those probes too)
    still_na <- is.na(imputed)
    n_knn_filled <- sum(!still_na)
    if (any(still_na)) {
      missing_probe_ids <- rownames(betas)[missing_mask][still_na]
      fallback <- ref_means[missing_probe_ids]
      fallback[is.na(fallback)] <- 0.5
      imputed[still_na] <- fallback
      diag_df$n_imputed_zeroshot[i] <- length(fallback)
    }

    betas[missing_mask, i] <- imputed
    diag_df$n_imputed_knn[i] <- n_knn_filled
    pb$tick()
  }

  list(betas = betas, sample_info = diag_df)
}


# ---------------------------------------------------------------------------
# Helper: per-sample reference mean fill
# ---------------------------------------------------------------------------

#' @keywords internal
impute_from_reference <- function(sample_vec, missing_mask, ref_means) {
  missing_probe_ids <- names(sample_vec)[missing_mask]
  fill <- if (!is.null(ref_means)) {
    fall <- ref_means[missing_probe_ids]
    fall[is.na(fall)] <- 0.5
    fall
  } else {
    rep(0.5, sum(missing_mask))
  }
  sample_vec[missing_mask] <- fill
  list(values = sample_vec, n_filled = sum(missing_mask))
}


# ---------------------------------------------------------------------------
# Helper: nearest-neighbour search within the input batch
# ---------------------------------------------------------------------------

#' Find k nearest neighbours of a target sample within the input batch
#'
#' Distance is computed only over probes that are non-missing in the target
#' AND in each candidate neighbour. Excludes the target column itself.
#'
#' @keywords internal
find_nearest_neighbors <- function(betas, target_idx, k,
                                     distance = c("pearson", "euclidean")) {
  distance <- match.arg(distance)
  n <- ncol(betas)
  if (n < 2L) return(NULL)

  target_vec <- betas[, target_idx]
  target_present <- !is.na(target_vec)
  if (sum(target_present) < 50L) return(NULL)

  # Subset rows once for speed; distance computed over target-present probes
  candidates <- setdiff(seq_len(n), target_idx)
  ref_block <- betas[target_present, candidates, drop = FALSE]
  tgt <- target_vec[target_present]

  d <- if (distance == "pearson") {
    cors <- suppressWarnings(stats::cor(tgt, ref_block,
                                          use = "pairwise.complete.obs"))
    cors <- as.numeric(cors)
    cors[!is.finite(cors)] <- -1
    1 - cors
  } else {
    diffs <- ref_block - tgt
    colMeans(diffs^2, na.rm = TRUE)
  }

  k_eff <- min(k, length(d))
  if (k_eff < 1L) return(NULL)
  ord <- order(d, decreasing = FALSE)[seq_len(k_eff)]
  list(idx = candidates[ord], dist = d[ord])
}


# ---------------------------------------------------------------------------
# Helper: NA-aware weighted row mean
# ---------------------------------------------------------------------------

#' @keywords internal
weighted_rowmean_with_na <- function(mat, weights) {
  if (anyNA(mat)) {
    not_na <- !is.na(mat)
    w_mat <- matrix(weights, nrow = nrow(mat), ncol = ncol(mat), byrow = TRUE)
    w_mat[!not_na] <- 0
    val_mat <- mat
    val_mat[!not_na] <- 0
    denom <- rowSums(w_mat)
    out <- rowSums(val_mat * w_mat) / denom
    out[denom == 0] <- NA_real_
    return(out)
  }
  drop(mat %*% weights)
}
