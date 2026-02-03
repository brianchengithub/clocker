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
#' @return Named vector of sex predictions (1 = Female, 0 = Male, 0.5 = Unknown)
#' @keywords internal
infer_sex_from_betas <- function(betas, platform = "EPIC", verbose = TRUE) {

  n_samples <- ncol(betas)
  sex_pred <- rep(0.5, n_samples)
  names(sex_pred) <- colnames(betas)

  # Get chromosome information from manifest
  probe_chr <- get_probe_chromosomes(platform, verbose = FALSE)

  if (is.null(probe_chr)) {
    if (verbose) message("    Could not load chromosome annotations")
    return(sex_pred)
  }

  common_probes <- intersect(names(probe_chr), rownames(betas))

  if (length(common_probes) == 0) {
    if (verbose) message("    No probes matched manifest")
    return(sex_pred)
  }

  x_probes <- common_probes[probe_chr[common_probes] %in% c("chrX", "X")]
  y_probes <- common_probes[probe_chr[common_probes] %in% c("chrY", "Y")]

  if (verbose) {
    message("    Found ", length(x_probes), " chrX probes, ", length(y_probes), " chrY probes")
  }

  if (length(x_probes) < 100 || length(y_probes) < 10) {
    if (verbose) message("    Insufficient sex chromosome probes for inference")
    return(sex_pred)
  }

  x_betas <- betas[x_probes, , drop = FALSE]
  y_betas <- betas[y_probes, , drop = FALSE]

  x_median <- apply(x_betas, 2, median, na.rm = TRUE)
  y_median <- apply(y_betas, 2, median, na.rm = TRUE)

  for (i in seq_len(n_samples)) {
    if (y_median[i] > 0.2) {
      sex_pred[i] <- 0       # Male
    } else if (y_median[i] < 0.1 && x_median[i] > 0.3) {
      sex_pred[i] <- 1       # Female
    } else {
      sex_pred[i] <- 0.5     # Ambiguous
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
