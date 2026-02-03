# ============================================================================
# Clock Computation Orchestrator
# Coordinates all clock calculations and assembles results
# ============================================================================


#' Check which clock packages are available
#' @param available_probes Character vector of available probe IDs
#' @param verbose Print progress
#' @return Named list of logical values indicating package availability
#' @keywords internal
check_clock_availability <- function(available_probes, verbose = TRUE) {

  availability <- list(
    sesame = requireNamespace("sesame", quietly = TRUE),
    dunedin_pace = requireNamespace("DunedinPACE", quietly = TRUE),
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
      if (pkg == "pc_clocks" && availability[[pkg]]) {
        status <- "available (via methylCIPHER)"
      }
      log_msg("  %s: %s", pkg, status, verbose = TRUE)
    }
  }

  return(availability)
}


#' Compute all clocks
#'
#' Main orchestrator that runs direct clock implementations, EpiDISH cell
#' deconvolution, sex inference, DunedinPACE, PC clocks, mitotic clocks,
#' and additional methylCIPHER clocks.
#'
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param pheno Optional phenotype data.frame
#' @param n_cores Number of CPU cores
#' @param verbose Print progress
#' @return Data frame with all clock results
#' @keywords internal
compute_all_clocks <- function(betas, pheno = NULL, n_cores, verbose = TRUE) {

  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)

  # ===== Direct weighted-sum clock implementations =====
  results <- compute_direct_clocks(betas, results, verbose)

  # ===== EpiDISH Cell Type Deconvolution =====
  results <- estimate_cell_composition(betas, results, verbose)

  # ===== Sex Inference =====
  results <- compute_sex_inference(betas, results, verbose)

  # ===== DunedinPACE =====
  results <- compute_dunedin_pace(betas, results, verbose)

  # ===== PC-Clocks via methylCIPHER =====
  results <- compute_pc_clocks(betas, results, pheno, verbose)

  # ===== EpiMitClocks (additional mitotic clocks) =====
  results <- compute_mitotic_clocks(betas, results, verbose)

  # ===== methylCIPHER additional clocks =====
  results <- compute_additional_clocks(betas, results, verbose)

  return(results)
}


# --- Individual clock computation helpers ---


#' Compute direct weighted-sum clocks
#' @keywords internal
compute_direct_clocks <- function(betas, results, verbose = TRUE) {

  if (verbose) message("  Computing clocks via direct implementations...")

  tryCatch({
    coeffs <- initialize_clock_coefficients()

    if (length(coeffs) > 0 && verbose) {
      message("    Loaded ", length(coeffs), " coefficient datasets")
    }

    mc_clocks <- list(
      "Horvath1"   = "Horvath1_CpGs",
      "Horvath2"   = "Horvath2_CpGs",
      "Hannum"     = "Hannum_CpGs",
      "PhenoAge"   = "PhenoAge_CpGs",
      "DNAmTL"     = "DNAmTL_CpGs",
      "Lin"        = "Lin_CpGs",
      "Zhang"      = "Zhang_10_CpG",
      "Zhang2019"  = "Zhang2019_CpGs",
      "Bocklandt"  = "Bocklandt_CpG",
      "Weidner"    = "Weidner_CpGs",
      "VidalBralo" = "VidalBralo_CpGs",
      "Garagnani"  = "Garagnani_CpG"
    )

    horvath_transform_clocks <- c("Horvath1", "Horvath2")
    computed_direct <- c()

    for (clock_name in names(mc_clocks)) {
      coef_name <- mc_clocks[[clock_name]]
      if (coef_name %in% names(coeffs)) {
        tryCatch({
          if (clock_name %in% horvath_transform_clocks) {
            result <- calc_weighted_sum_clock(betas, coeffs[[coef_name]],
                                              transform_func = horvath_age_transform)
          } else {
            result <- calc_weighted_sum_clock(betas, coeffs[[coef_name]])
          }
          if (!is.null(result) && length(result) == ncol(betas)) {
            results[[clock_name]] <- as.numeric(result)
            computed_direct <- c(computed_direct, clock_name)
          } else if (verbose && is.null(result)) {
            message("    ", clock_name, ": no matching probes found")
          }
        }, error = function(e) {
          if (verbose) message("    ", clock_name, " error: ", e$message)
        })
      }
    }

    # epiTOC2 direct calculation
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

  return(results)
}


#' Compute sex inference and add to results
#' @keywords internal
compute_sex_inference <- function(betas, results, verbose = TRUE) {

  tryCatch({
    log_msg("  Inferring sex...", verbose = verbose)

    n_probes <- nrow(betas)
    platform <- if (n_probes > 900000) "EPICv2"
                else if (n_probes > 800000) "EPIC"
                else if (n_probes > 400000) "HM450"
                else "EPIC"

    sex_numeric <- infer_sex_from_betas(betas, platform = platform, verbose = verbose)

    if (!is.null(sex_numeric) && length(sex_numeric) == ncol(betas)) {
      results$InferredSex_Numeric <- sex_numeric
      results$InferredSex <- ifelse(sex_numeric == 1, "F",
                                    ifelse(sex_numeric == 0, "M", "U"))
    }
  }, error = function(e) {
    if (verbose) message("    Sex inference failed: ", e$message)
  })

  return(results)
}


#' Compute DunedinPACE
#' @keywords internal
compute_dunedin_pace <- function(betas, results, verbose = TRUE) {

  if (!requireNamespace("DunedinPACE", quietly = TRUE)) return(results)

  tryCatch({
    log_msg("  Calculating DunedinPACE...", verbose = verbose)

    # DunedinPACE expects probes as rows, samples as columns
    # Ensure it's a proper matrix (critical for single-sample case)
    betas_dp <- betas
    if (!is.matrix(betas_dp)) {
      betas_dp <- as.matrix(betas_dp)
    }

    # Force 2D even for single sample
    if (is.null(dim(betas_dp))) {
      probe_names <- names(betas_dp)
      betas_dp <- matrix(betas_dp, ncol = 1,
                         dimnames = list(probe_names, colnames(betas)[1]))
    }

    pace <- DunedinPACE::PACEProjector(betas_dp)

    if (!is.null(pace)) {
      pace_val <- NULL
      if (is.list(pace) && "DunedinPACE" %in% names(pace)) {
        pace_val <- pace$DunedinPACE
      } else if (is.data.frame(pace) && "DunedinPACE" %in% colnames(pace)) {
        pace_val <- pace$DunedinPACE
      } else if (is.numeric(pace)) {
        pace_val <- pace
      }

      if (!is.null(pace_val)) {
        # Handle single-value result
        if (length(pace_val) == 1 && ncol(betas) == 1) {
          results$DunedinPACE <- pace_val
        } else if (length(pace_val) == ncol(betas)) {
          results$DunedinPACE <- pace_val
        }
      }
    }
  }, error = function(e) {
    if (verbose) message("    DunedinPACE failed: ", e$message)
  })

  return(results)
}


#' Compute PC clocks via methylCIPHER
#' @keywords internal
compute_pc_clocks <- function(betas, results, pheno = NULL, verbose = TRUE) {

  if (!requireNamespace("methylCIPHER", quietly = TRUE)) return(results)

  tryCatch({
    if (verbose) message("  Computing PC clocks...")

    # Check for and download PC clocks data file if needed
    pc_data_path <- NULL
    cache_dir <- get_manifest_cache_dir()
    pc_data_file <- file.path(cache_dir, "PCClocks_data.qs2")

    if (file.exists(pc_data_file)) {
      file_size <- file.info(pc_data_file)$size
      if (file_size > 1e9) {
        pc_data_path <- pc_data_file
        if (verbose) message("    Using cached PC clocks data (", round(file_size/1e9, 2), " GB)")
      } else {
        if (verbose) message("    Cached PC clocks data appears incomplete, re-downloading...")
        unlink(pc_data_file)
      }
    }

    if (is.null(pc_data_path)) {
      if (verbose) message("    Downloading PC clocks data from Zenodo (~2GB, this may take several minutes)...")
      pc_url <- "https://zenodo.org/records/17162604/files/PCClocks_data.qs2?download=1"

      old_timeout <- getOption("timeout")
      options(timeout = 1800)

      tryCatch({
        download.file(pc_url, pc_data_file, mode = "wb", quiet = !verbose)
      }, error = function(e) {
        if (verbose) message("    Failed to download PC clocks data: ", e$message)
      }, warning = function(w) {
        if (verbose) message("    Download warning: ", w$message)
      })

      options(timeout = old_timeout)

      if (file.exists(pc_data_file)) {
        file_size <- file.info(pc_data_file)$size
        if (file_size > 1e9) {
          pc_data_path <- pc_data_file
          if (verbose) message("    PC clocks data downloaded and cached (", round(file_size/1e9, 2), " GB)")
        } else {
          if (verbose) message("    Download incomplete, removing partial file")
          unlink(pc_data_file)
        }
      }
    }

    # Load PC clocks CpG data to global env
    tryCatch({
      data("PCClocks_CpGs", package = "methylCIPHER", envir = .GlobalEnv)
    }, error = function(e) NULL, warning = function(w) NULL)

    betas_t <- t(betas)
    sample_ids <- rownames(betas_t)

    if (!exists("calcPCClocks", envir = asNamespace("methylCIPHER"))) {
      if (verbose) message("    calcPCClocks function not found in methylCIPHER")
      return(results)
    }

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

    pheno_df <- data.frame(
      Sample_ID = sample_ids,
      Age = age_values,
      Female = female_values,
      stringsAsFactors = FALSE
    )
    rownames(pheno_df) <- sample_ids

    if (verbose) {
      message("    PC clocks pheno: Age range = ",
              round(min(pheno_df$Age, na.rm = TRUE), 1), " - ",
              round(max(pheno_df$Age, na.rm = TRUE), 1),
              ", Female: ", sum(pheno_df$Female == 1), " F, ",
              sum(pheno_df$Female == 0), " M")
    }

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
          message("    PC result cols: ", paste(head(colnames(pc_result), 10), collapse = ", "))
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
      }
    } else {
      if (verbose) message("    PC clocks: No result returned")
    }

  }, error = function(e) {
    if (verbose) message("    PC clocks error: ", e$message)
  })

  return(results)
}


#' Compute mitotic clocks via EpiMitClocks
#' @keywords internal
compute_mitotic_clocks <- function(betas, results, verbose = TRUE) {

  epi_pkg <- NULL
  for (pkg_name in c("EpiMitClocks", "epiTOC2")) {
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      epi_pkg <- pkg_name
      break
    }
  }

  if (is.null(epi_pkg) || "epiTOC2_TNSC" %in% colnames(results)) return(results)

  tryCatch({
    if (verbose) message("  Computing mitotic clocks via EpiMitClocks...")

    epi_data_items <- c("dataETOC3", "cugpmitclockCpG", "epiTOCcpgs3", "estETOC3",
                        "EpiCMITcpgs", "Replitali")
    for (d in epi_data_items) {
      tryCatch({
        data(list = d, package = epi_pkg, envir = .GlobalEnv)
      }, error = function(e) NULL, warning = function(w) NULL)
    }

    if (exists("EpiMitClocks", envir = asNamespace(epi_pkg))) {
      tryCatch({
        epi_results <- get("EpiMitClocks", envir = asNamespace(epi_pkg))(
          data.m = betas, ages.v = NULL
        )

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

  return(results)
}


#' Compute additional methylCIPHER clocks
#' @keywords internal
compute_additional_clocks <- function(betas, results, verbose = TRUE) {

  if (!requireNamespace("methylCIPHER", quietly = TRUE)) return(results)

  tryCatch({
    if (verbose) message("  Computing additional clocks via methylCIPHER functions...")

    betas_t <- t(betas)

    additional_clocks <- c("AdaptAge", "CausAge", "DamAge", "SystemsAge")

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

            if (!is.null(result) && is.numeric(result)) {
              # Handle single-sample: result might be length 1
              if (length(result) == ncol(betas) ||
                  (length(result) == 1 && ncol(betas) == 1)) {
                results[[clock]] <- as.numeric(result)
                added_clocks <- c(added_clocks, clock)
              }
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

  return(results)
}
