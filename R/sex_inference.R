# ============================================================================
# Sex Inference (exact methylQC algorithm)
#
# This file implements the sex caller from the methylQC R package
# (https://github.com/brianchengithub/methylQC, MIT license, by Brian Cheng).
# The algorithm is reproduced verbatim from methylQC's
# `compute_sex_intensities_single()` and `plot_sex_check_optimal()` functions
# so that clocker's sex calls are identical to a methylQC run on the
# same input.
#
# Algorithm:
#
#   Stage 1 -- per-sample sex signals
#     Two curated probe sets from sesameData::sesameDataGet("EPIC.probeInfo"):
#       * chrX_xlinked  (3,433 probes): non-PAR X-inactivated chrX probes
#       * chrY_clean    (   314 probes): non-PAR chrY probes minus cross-hybs
#     Per sample:
#       chrX_signal = median over chrX_xlinked probes of total signal intensity
#                     ( MG + MR + UG + UR, with NA -> 0 ) when SigDFs are
#                     available. When only a beta matrix is supplied,
#                     chrX_signal = median beta over the same probes.
#       chrY_signal: analogous using chrY_clean.
#     Both require >= 10 matching probes; otherwise NA.
#
#   Stage 2 -- data-driven threshold optimization
#     Candidate thresholds: tau in quantile(chrY, seq(0.15, 0.85, 0.01))
#     For each tau:
#       females = samples with chrY <= tau   (require n >= 3)
#       males   = samples with chrY >  tau   (require n >= 3)
#       fit lm(chrY ~ chrX) in each cluster
#       cost(tau) = sum |residuals_F| + sum |residuals_M|
#     Best tau = argmin cost (default: median(chrY) if no valid tau).
#
#   Stage 3 -- final cluster regressions and orthogonal-distance bands
#     Refit lm(chrY ~ chrX) on the final F and M clusters.
#     Orthogonal distance from a point (x,y) to a fitted line (slope a,
#     intercept c):  d(x,y) = | a*x - y + c | / sqrt(a^2 + 1).
#     sd_F = sd(d_F at F samples); sd_M = sd(d_M at M samples)
#     Band thresholds: thr_F = 5 * sd_F ; thr_M = 5 * sd_M
#
#   Stage 4 -- classification
#     d_F <= thr_F  &  d_M >  thr_M  -> "F"
#     d_M <= thr_M  &  d_F >  thr_F  -> "M"
#     d_F <= thr_F  &  d_M <= thr_M  -> tiebreak by chrY threshold
#                                        (chrY <= tau -> F else M)
#     otherwise                       -> "Unclear"
#
# Small-batch behaviour:
#   methylQC's algorithm requires >= 6 input samples. clocker may be
#   run on a single sample or a small batch where the data-driven cluster
#   fit is impossible. In that case we fall back to the embedded
#   reference-cluster regressions in DEFAULT_SEX_REFERENCE (calibrated on a
#   blood-derived training set; user-tunable via set_sex_reference()) and
#   apply the same Stage 3 / Stage 4 logic with the reference parameters.
# ============================================================================


# ---------------------------------------------------------------------------
# Reference cluster parameters for small-batch fallback
# ---------------------------------------------------------------------------

DEFAULT_SEX_REFERENCE <- list(
  intensity = list(
    threshold = 6500,    # chrY total-intensity cut between clusters
    male   = list(slope = -0.10, intercept = 7800, sigma = 350),
    female = list(slope = -0.02, intercept = 1100, sigma = 300)
  ),
  beta = list(
    threshold = 0.10,
    male   = list(slope = -0.10, intercept = 0.350, sigma = 0.030),
    female = list(slope = -0.02, intercept = 0.040, sigma = 0.020)
  ),
  notes = paste(
    "Defaults calibrated on whole-blood EPIC reference cohorts.",
    "Override with set_sex_reference() for atypical tissues."
  )
)


#' Set custom reference parameters for the small-batch sex caller
#' @param scale "intensity" (IDAT input) or "beta" (matrix input)
#' @param threshold chrY threshold separating F (below) from M (above)
#' @param male,female Lists with `slope`, `intercept`, `sigma`
#' @export
set_sex_reference <- function(scale = c("intensity", "beta"),
                                threshold = NULL, male = NULL, female = NULL) {
  scale <- match.arg(scale)
  cur <- .qc_env$sex_reference %||% DEFAULT_SEX_REFERENCE
  if (!is.null(threshold)) cur[[scale]]$threshold <- threshold
  if (!is.null(male))      cur[[scale]]$male      <- male
  if (!is.null(female))    cur[[scale]]$female    <- female
  .qc_env$sex_reference <- cur
  invisible(NULL)
}


#' @keywords internal
get_sex_reference <- function() {
  .qc_env$sex_reference %||% DEFAULT_SEX_REFERENCE
}


# ---------------------------------------------------------------------------
# Stage 1 -- per-sample sex signals
# ---------------------------------------------------------------------------

#' Compute (chrX_signal, chrY_signal) for a single SeSAMe SigDF
#'
#' Reproduces methylQC::compute_sex_intensities_single() exactly.
#'
#' @keywords internal
compute_sex_intensities_single <- function(sdf) {
  pid <- sdf$Probe_ID
  ti <- rowSums(cbind(
          ifelse(is.na(sdf$MG), 0, sdf$MG),
          ifelse(is.na(sdf$MR), 0, sdf$MR),
          ifelse(is.na(sdf$UG), 0, sdf$UG),
          ifelse(is.na(sdf$UR), 0, sdf$UR)
        ), na.rm = TRUE)
  xi <- which(pid %in% .sesame_chrX_xlinked)
  yi <- which(pid %in% .sesame_chrY_clean)
  list(
    chrX = if (length(xi) >= 10) stats::median(ti[xi], na.rm = TRUE) else NA_real_,
    chrY = if (length(yi) >= 10) stats::median(ti[yi], na.rm = TRUE) else NA_real_
  )
}


#' Compute sex signals from a list of SigDFs
#' @keywords internal
compute_sex_signals_from_sdfs <- function(sdfs) {
  if (is.null(sdfs) || length(sdfs) == 0L) return(NULL)
  rows <- lapply(seq_along(sdfs), function(i) {
    si <- compute_sex_intensities_single(sdfs[[i]])
    data.frame(sample_id = names(sdfs)[i] %||% paste0("Sample_", i),
                chrX = si$chrX, chrY = si$chrY,
                stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}


#' Compute sex signals from a beta matrix (proxy when SigDFs not available)
#'
#' For beta-only input we substitute median beta over the same curated
#' probe sets in place of total intensity. The downstream classification
#' algorithm is identical; only the per-sample summary statistic differs.
#'
#' @keywords internal
compute_sex_signals_from_betas <- function(betas) {
  pid <- rownames(betas)
  xi <- which(pid %in% .sesame_chrX_xlinked)
  yi <- which(pid %in% .sesame_chrY_clean)

  if (length(xi) < 10 || length(yi) < 10) {
    return(data.frame(
      sample_id = colnames(betas),
      chrX = NA_real_, chrY = NA_real_,
      stringsAsFactors = FALSE))
  }

  if (requireNamespace("matrixStats", quietly = TRUE)) {
    chrX_vals <- matrixStats::colMedians(betas[xi, , drop = FALSE], na.rm = TRUE)
    chrY_vals <- matrixStats::colMedians(betas[yi, , drop = FALSE], na.rm = TRUE)
  } else {
    chrX_vals <- apply(betas[xi, , drop = FALSE], 2, stats::median, na.rm = TRUE)
    chrY_vals <- apply(betas[yi, , drop = FALSE], 2, stats::median, na.rm = TRUE)
  }
  data.frame(
    sample_id = colnames(betas),
    chrX = as.numeric(chrX_vals),
    chrY = as.numeric(chrY_vals),
    stringsAsFactors = FALSE)
}


# ---------------------------------------------------------------------------
# Stage 2 + 3 + 4 -- methylQC clustering algorithm
# Reproduces methylQC::plot_sex_check_optimal() exactly (sans plotting).
# ---------------------------------------------------------------------------

#' Orthogonal distance from points to a fitted line
#' @keywords internal
ortho_dist <- function(x, y, fit) {
  if (is.null(fit)) return(rep(Inf, length(x)))
  a <- stats::coef(fit)[2]; c <- stats::coef(fit)[1]
  abs(a * x - y + c) / sqrt(a^2 + 1)
}


#' Run the methylQC sex-cluster algorithm
#'
#' Verbatim implementation of methylQC's plot_sex_check_optimal() classifier.
#' Requires at least 6 samples; smaller batches return NULL so the caller
#' can fall back to the reference-based classifier.
#'
#' @keywords internal
classify_methylqc_clusters <- function(sex_df, verbose = TRUE) {

  if (nrow(sex_df) < 6L) return(NULL)

  candidates <- stats::quantile(sex_df$chrY, probs = seq(0.15, 0.85, by = 0.01),
                                  na.rm = TRUE)
  best_thresh <- stats::median(sex_df$chrY)
  best_cost <- Inf
  for (thr in candidates) {
    fi <- which(sex_df$chrY <= thr)
    mi <- which(sex_df$chrY >  thr)
    if (length(fi) < 3 || length(mi) < 3) next
    ff <- tryCatch(stats::lm(chrY ~ chrX, data = sex_df[fi, ]),
                    error = function(e) NULL)
    fm <- tryCatch(stats::lm(chrY ~ chrX, data = sex_df[mi, ]),
                    error = function(e) NULL)
    if (is.null(ff) || is.null(fm)) next
    cost <- sum(abs(stats::residuals(ff))) + sum(abs(stats::residuals(fm)))
    if (cost < best_cost) { best_cost <- cost; best_thresh <- thr }
  }

  f_idx <- which(sex_df$chrY <= best_thresh)
  m_idx <- which(sex_df$chrY >  best_thresh)
  fit_f <- if (length(f_idx) >= 3) tryCatch(
    stats::lm(chrY ~ chrX, data = sex_df[f_idx, ]),
    error = function(e) NULL) else NULL
  fit_m <- if (length(m_idx) >= 3) tryCatch(
    stats::lm(chrY ~ chrX, data = sex_df[m_idx, ]),
    error = function(e) NULL) else NULL

  if (is.null(fit_f) && is.null(fit_m)) return(NULL)

  d_f <- ortho_dist(sex_df$chrX, sex_df$chrY, fit_f)
  d_m <- ortho_dist(sex_df$chrX, sex_df$chrY, fit_m)
  sd_f <- if (length(f_idx) >= 3) stats::sd(d_f[f_idx]) else Inf
  sd_m <- if (length(m_idx) >= 3) stats::sd(d_m[m_idx]) else Inf

  band_sd <- 5.0
  thresh_f <- band_sd * sd_f
  thresh_m <- band_sd * sd_m

  sex_df$inferred_sex_intensity <- ifelse(
    d_f <= thresh_f & d_m >  thresh_m, "F",
    ifelse(d_m <= thresh_m & d_f >  thresh_f, "M",
    ifelse(d_f <= thresh_f & d_m <= thresh_m,
           ifelse(sex_df$chrY <= best_thresh, "F", "M"), "Unclear")))

  sex_df$flag <- ifelse(
    d_f <= thresh_f & d_m >  thresh_m, "F_band_only",
    ifelse(d_m <= thresh_m & d_f >  thresh_f, "M_band_only",
    ifelse(d_f <= thresh_f & d_m <= thresh_m,
           "both_bands_chrY_tiebreak", "outside_both_bands")))

  attr(sex_df, "best_threshold") <- as.numeric(best_thresh)
  attr(sex_df, "fit_f") <- fit_f
  attr(sex_df, "fit_m") <- fit_m
  attr(sex_df, "sd_f")  <- sd_f
  attr(sex_df, "sd_m")  <- sd_m
  attr(sex_df, "n_f")   <- length(f_idx)
  attr(sex_df, "n_m")   <- length(m_idx)
  attr(sex_df, "cluster_cost")    <- best_cost
  attr(sex_df, "used_reference")  <- FALSE

  if (verbose) {
    log_msg(paste("    Sex clustering: threshold=%.3f, F=%d, M=%d,",
                    "unclear=%d (chrY range %.3f-%.3f)"),
              best_thresh,
              attr(sex_df, "n_f"), attr(sex_df, "n_m"),
              sum(sex_df$inferred_sex_intensity == "Unclear"),
              min(sex_df$chrY, na.rm = TRUE),
              max(sex_df$chrY, na.rm = TRUE),
              verbose = TRUE)
  }
  sex_df
}


#' Reference-based fallback (n < 6 or pathological cluster fit)
#'
#' Uses embedded reference parameters (DEFAULT_SEX_REFERENCE) for the
#' threshold and per-cluster regressions, then runs the same Stage 3/4
#' logic so the classification semantics match a real methylQC fit.
#'
#' @keywords internal
classify_with_reference <- function(sex_df, scale = c("intensity", "beta"),
                                      verbose = TRUE) {
  scale <- match.arg(scale)
  ref <- get_sex_reference()[[scale]]

  make_fit <- function(slope, intercept) {
    fit <- list()
    fit$coefficients <- c("(Intercept)" = intercept, "chrX" = slope)
    class(fit) <- "lm"
    fit
  }
  fit_f <- make_fit(ref$female$slope, ref$female$intercept)
  fit_m <- make_fit(ref$male$slope,   ref$male$intercept)

  d_f <- ortho_dist(sex_df$chrX, sex_df$chrY, fit_f)
  d_m <- ortho_dist(sex_df$chrX, sex_df$chrY, fit_m)
  thresh_f <- 5.0 * ref$female$sigma
  thresh_m <- 5.0 * ref$male$sigma

  sex_df$inferred_sex_intensity <- ifelse(
    d_f <= thresh_f & d_m >  thresh_m, "F",
    ifelse(d_m <= thresh_m & d_f >  thresh_f, "M",
    ifelse(d_f <= thresh_f & d_m <= thresh_m,
           ifelse(sex_df$chrY <= ref$threshold, "F", "M"), "Unclear")))

  sex_df$flag <- ifelse(
    d_f <= thresh_f & d_m >  thresh_m, "F_band_only_ref",
    ifelse(d_m <= thresh_m & d_f >  thresh_f, "M_band_only_ref",
    ifelse(d_f <= thresh_f & d_m <= thresh_m,
           "both_bands_chrY_tiebreak_ref", "outside_both_bands_ref")))

  attr(sex_df, "best_threshold") <- ref$threshold
  attr(sex_df, "fit_f") <- fit_f
  attr(sex_df, "fit_m") <- fit_m
  attr(sex_df, "sd_f")  <- ref$female$sigma
  attr(sex_df, "sd_m")  <- ref$male$sigma
  attr(sex_df, "scale") <- scale
  attr(sex_df, "used_reference") <- TRUE

  if (verbose) {
    log_msg(paste("    Sex clustering (reference fallback, n < 6):",
                    "scale=%s, threshold=%.3f"),
              scale, ref$threshold, verbose = TRUE)
  }
  sex_df
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

#' Infer sex from DNA methylation using the methylQC algorithm
#'
#' Two input modes:
#'  * `sdfs` not NULL : intensities are computed exactly as in methylQC
#'                      using sum(MG, MR, UG, UR). Preferred mode; used
#'                      automatically when clocker() reads from IDATs.
#'  * `sdfs` NULL     : a beta-value proxy is used (median beta over the
#'                      same curated probe sets). The classification
#'                      algorithm is otherwise identical.
#'
#' Inputs may also come pre-cached in `.qc_env$sex_signals` (set by
#' input_processing.R when SigDFs are seen but not retained in memory).
#'
#' @keywords internal
infer_sex <- function(betas, sdfs = NULL, verbose = TRUE) {

  n_samples <- ncol(betas)
  empty <- list(
    chrX_signal = setNames(rep(NA_real_, n_samples), colnames(betas)),
    chrY_signal = setNames(rep(NA_real_, n_samples), colnames(betas)),
    sex         = setNames(rep("Unclear", n_samples), colnames(betas)),
    flag        = setNames(rep("uncalled", n_samples), colnames(betas)),
    scale       = NA_character_,
    details     = NULL
  )

  sex_signals <- .qc_env$sex_signals
  scale <- "beta"

  if (!is.null(sex_signals) && nrow(sex_signals) == n_samples) {
    scale <- attr(sex_signals, "scale") %||% "intensity"
    if (verbose) {
      log_msg("    Sex inference: using cached %s signals (n=%d)",
              scale, nrow(sex_signals), verbose = TRUE)
    }
  } else if (!is.null(sdfs)) {
    if (verbose) {
      log_msg("    Sex inference: total-intensity signals from SigDFs (n=%d)",
              length(sdfs), verbose = TRUE)
    }
    sex_signals <- compute_sex_signals_from_sdfs(sdfs)
    attr(sex_signals, "scale") <- "intensity"
    scale <- "intensity"
  } else {
    if (verbose) {
      log_msg("    Sex inference: no SigDFs available -- using beta-value proxy",
              verbose = TRUE)
    }
    sex_signals <- compute_sex_signals_from_betas(betas)
    attr(sex_signals, "scale") <- "beta"
    scale <- "beta"
  }

  if (!identical(sex_signals$sample_id, colnames(betas))) {
    sex_signals <- sex_signals[match(colnames(betas), sex_signals$sample_id), ,
                                drop = FALSE]
    sex_signals$sample_id[is.na(sex_signals$sample_id)] <-
      colnames(betas)[is.na(sex_signals$sample_id)]
  }

  complete <- !is.na(sex_signals$chrX) & !is.na(sex_signals$chrY)
  if (!any(complete)) {
    if (verbose) message("    Sex inference: no complete (chrX, chrY) pairs available")
    empty$chrX_signal <- setNames(sex_signals$chrX, sex_signals$sample_id)
    empty$chrY_signal <- setNames(sex_signals$chrY, sex_signals$sample_id)
    empty$scale <- scale
    return(empty)
  }

  classify_input <- sex_signals[complete, , drop = FALSE]

  classified <- classify_methylqc_clusters(classify_input, verbose = verbose)
  if (is.null(classified)) {
    classified <- classify_with_reference(classify_input, scale = scale,
                                            verbose = verbose)
  }

  out_sex  <- setNames(rep("Unclear",  nrow(sex_signals)), sex_signals$sample_id)
  out_flag <- setNames(rep("missing_signal", nrow(sex_signals)), sex_signals$sample_id)
  out_sex[classified$sample_id]  <- classified$inferred_sex_intensity
  out_flag[classified$sample_id] <- classified$flag

  out_sex  <- out_sex[colnames(betas)]
  out_flag <- out_flag[colnames(betas)]

  if (verbose) {
    cnt <- table(factor(out_sex, levels = c("F", "M", "Unclear")))
    log_msg("    Inferred sex: F=%d, M=%d, Unclear=%d",
            cnt["F"], cnt["M"], cnt["Unclear"], verbose = TRUE)
  }

  list(
    chrX_signal = setNames(sex_signals$chrX, sex_signals$sample_id)[colnames(betas)],
    chrY_signal = setNames(sex_signals$chrY, sex_signals$sample_id)[colnames(betas)],
    sex         = out_sex,
    flag        = out_flag,
    scale       = scale,
    details     = classified
  )
}


#' Cache pre-computed sex signals (called from input_processing.R when
#' SigDFs are processed in streaming fashion so we don't have to retain
#' SigDFs in memory).
#' @keywords internal
cache_sex_signals <- function(sex_signal_df, scale = c("intensity", "beta")) {
  scale <- match.arg(scale)
  attr(sex_signal_df, "scale") <- scale
  .qc_env$sex_signals <- sex_signal_df
  invisible(NULL)
}
