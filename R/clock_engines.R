# ============================================================================
# Clock Calculation Engines
# Low-level functions: coefficient loading, weighted sum, direct calculations
# ============================================================================


#' Initialize clock coefficients from installed packages
#'
#' Loads coefficient datasets from methylCIPHER and EpiMitClocks into memory.
#' Uses isolated environments to avoid polluting the global namespace.
#'
#' @return Named list of coefficient data frames
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
      "dataETOC3"       = "dataETOC3.l",
      "estETOC3"        = "estETOC3.m",
      "epiTOCcpgs3"     = "epiTOCcpgs3.v",
      "cugpmitclockCpG" = "cugpmitclockCpG.v",
      "EpiCMITcpgs"     = "epiCMIT.df",
      "Replitali"       = c("replitali.coe", "replitali.cpg.v")
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
#'
#' Converts raw clock values to DNAm age in years using the Horvath
#' piecewise transformation.
#'
#' @param x Raw clock value (weighted sum)
#' @param adult_age Adult age constant (default 20 for Horvath clocks)
#' @return DNAm age in years
#' @keywords internal
horvath_age_transform <- function(x, adult_age = 20) {
  ifelse(x < 0,
         (adult_age + 1) * exp(x) - 1,
         (adult_age + 1) * x + adult_age)
}


#' Generic weighted sum clock calculator
#'
#' Computes a clock value as intercept + sum(beta * weight) with automatic
#' detection of CpG and weight columns in coefficient data frames.
#' Handles edge cases like single-CpG clocks, named vectors, and
#' non-standard column names.
#'
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coef_df Data frame, matrix, or named vector with CpG/weight info
#' @param transform_func Optional transformation function to apply to results
#' @return Named vector of clock values, or NULL on failure
#' @keywords internal
calc_weighted_sum_clock <- function(betas, coef_df, transform_func = NULL) {

  # Handle non-data.frame inputs (named vectors, single values, etc.)
  if (is.null(coef_df)) return(NULL)

  if (is.numeric(coef_df) && !is.matrix(coef_df) && !is.data.frame(coef_df)) {
    # Named numeric vector: names are CpGs, values are weights
    if (!is.null(names(coef_df))) {
      coef_df <- data.frame(CpG = names(coef_df), Coefficient = as.numeric(coef_df),
                            stringsAsFactors = FALSE)
    } else {
      return(NULL)
    }
  }

  if (is.matrix(coef_df)) {
    coef_df <- as.data.frame(coef_df, stringsAsFactors = FALSE)
    # If matrix had rownames but no CpG column, use rownames
    if (!is.null(rownames(coef_df)) && !any(sapply(coef_df, is.character))) {
      coef_df$CpG <- rownames(coef_df)
    }
  }

  if (!is.data.frame(coef_df) || ncol(coef_df) == 0) return(NULL)
  if (nrow(coef_df) == 0) return(NULL)

  # Auto-detect CpG column
  cpg_candidates <- c("CpG", "CpGmarker", "probe", "Probe", "cpg", "ID",
                       "ProbeID", "probe_id", "CpG_ID", "Marker")
  cpg_col <- intersect(cpg_candidates, colnames(coef_df))[1]
  if (is.na(cpg_col)) {
    # Try first character column
    char_cols <- which(sapply(coef_df, is.character) | sapply(coef_df, is.factor))
    if (length(char_cols) > 0) {
      cpg_col <- colnames(coef_df)[char_cols[1]]
    } else if (!is.null(rownames(coef_df)) &&
               any(grepl("^cg|^ch", rownames(coef_df), ignore.case = TRUE))) {
      # Use rownames as CpG IDs
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
    # Find first numeric column that's not the CpG column
    numeric_cols <- which(sapply(coef_df, is.numeric))
    numeric_cols <- setdiff(names(numeric_cols), cpg_col)
    if (length(numeric_cols) > 0) {
      weight_col <- numeric_cols[1]
    }
  }

  if (is.na(weight_col) || is.null(weight_col)) return(NULL)

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

  if (length(cpgs) == 0) return(NULL)

  matched_idx <- match(cpgs, rownames(betas))
  valid <- !is.na(matched_idx)

  if (sum(valid) == 0) return(NULL)

  betas_subset <- betas[matched_idx[valid], , drop = FALSE]
  weights_valid <- weights[valid]

  clock_values <- intercept + colSums(betas_subset * weights_valid, na.rm = TRUE)

  if (!is.null(transform_func)) {
    clock_values <- transform_func(clock_values)
  }

  return(clock_values)
}


#' Calculate epiTOC2 mitotic clock directly
#'
#' Computes the total number of stem cell divisions (TNSC) using the
#' epiTOC2 algorithm.
#'
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named vector of TNSC values, or NULL on failure
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
#'
#' Placeholder for direct PC clock calculation. Currently returns NULL
#' because PC clocks require the full training data from Zenodo.
#' The methylCIPHER::calcPCClocks function is used instead.
#'
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param coeffs Pre-loaded coefficients list
#' @return Named list of PC clock values, or NULL
#' @keywords internal
calc_pcclocks_direct <- function(betas, coeffs) {

  if (!"PCClocks_CpGs" %in% names(coeffs)) return(NULL)

  # PC clocks require the full PC training data (~2GB)
  # Use methylCIPHER::calcPCClocks with the downloaded data file instead
  return(NULL)
}
