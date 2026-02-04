# ============================================================================
# Sex Inference
# Infer biological sex from X/Y chromosome methylation patterns
# ============================================================================


#' Infer sex from DNA methylation data using X/Y chromosome probes
#'
#' Uses median methylation on X and Y chromosomes to predict biological sex.
#' Females (XX) show higher X methylation due to X-inactivation and low Y signal.
#' Males (XY) show lower X methylation and higher Y signal.
#'
#' @param betas Beta value matrix (CpGs as rows, samples as columns)
#' @param platform Platform name (EPIC, HM450, EPICv2, etc.)
#' @param verbose Print progress
#' @return Named list with: x_median, y_median (numeric vectors), sex (character vector)
#' @keywords internal
infer_sex_from_betas <- function(betas, platform = "EPIC", verbose = TRUE) {

  n_samples <- ncol(betas)
  result <- list(
    x_median = rep(NA_real_, n_samples),
    y_median = rep(NA_real_, n_samples),
    sex = rep("U", n_samples)
  )
  names(result$x_median) <- colnames(betas)
  names(result$y_median) <- colnames(betas)
  names(result$sex) <- colnames(betas)

  # Get chromosome information from manifest
  # After EPICv2 normalization, probes are base CpG names (e.g., cg00000029).
  # EPIC/HM450 manifests use base names → best match. EPICv2 manifest uses
  # suffixed names → only ~1500 non-CpG probes match, giving misleading results.
  # So try EPIC first for EPICv2 data.
  probe_chr <- NULL
  platforms_to_try <- platform
  if (platform %in% c("EPICv2", "EPICv2/EPIC+")) {
    platforms_to_try <- c("EPIC", "HM450", "EPICv2")
  }

  for (plt in platforms_to_try) {
    probe_chr <- get_probe_chromosomes(plt, verbose = FALSE)
    if (!is.null(probe_chr)) {
      common <- length(intersect(names(probe_chr), rownames(betas)))
      if (verbose) message("    Trying ", plt, " manifest: ", common, " common probes")
      if (common > 10000) break  # Need substantial overlap
      probe_chr <- NULL  # Not enough overlap, try next
    }
  }

  if (is.null(probe_chr)) {
    if (verbose) message("    Could not load chromosome annotations")
    return(result)
  }

  common_probes <- intersect(names(probe_chr), rownames(betas))

  if (length(common_probes) == 0) {
    if (verbose) message("    No probes matched manifest")
    return(result)
  }

  x_probes <- common_probes[probe_chr[common_probes] %in% c("chrX", "X")]
  y_probes <- common_probes[probe_chr[common_probes] %in% c("chrY", "Y")]

  if (verbose) {
    message("    Found ", length(x_probes), " chrX probes, ", length(y_probes), " chrY probes")
  }

  if (length(x_probes) < 100 || length(y_probes) < 10) {
    if (verbose) message("    Insufficient sex chromosome probes for inference")
    return(result)
  }

  x_betas <- betas[x_probes, , drop = FALSE]
  y_betas <- betas[y_probes, , drop = FALSE]

  result$x_median <- apply(x_betas, 2, median, na.rm = TRUE)
  result$y_median <- apply(y_betas, 2, median, na.rm = TRUE)

  for (i in seq_len(n_samples)) {
    if (result$y_median[i] > 0.2) {
      result$sex[i] <- "M"
    } else if (result$y_median[i] < 0.1 && result$x_median[i] > 0.3) {
      result$sex[i] <- "F"
    } else {
      result$sex[i] <- "U"
    }
  }

  if (verbose) {
    n_male <- sum(result$sex == "M")
    n_female <- sum(result$sex == "F")
    n_unknown <- sum(result$sex == "U")
    message("    Sex inference: ", n_female, " Female, ", n_male, " Male, ", n_unknown, " Unknown")
  }

  return(result)
}
