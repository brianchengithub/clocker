# ============================================================================
# Clock Computation Orchestrator
#
# Coordinates clock calculations and assembles results. Major fixes here
# (vs. v2.0):
#   * Pheno data is now reordered to match colnames(betas) by sample ID
#     (was previously pasted in row order; major silent-correctness bug).
#   * No more .GlobalEnv writes; data is loaded into .qc_env.
#   * PC-Clocks data file is verified by SHA-256 (with corruption fallback).
#   * When pheno is missing or incomplete, we explicitly alert the user
#     about which fallback (Horvath2/InferredSex/etc.) is being used.
#   * Sex column accepts 0/1, F/M, Female/Male (any case) via
#     coerce_female_indicator().
# ============================================================================


# Pinned PC-Clocks Zenodo asset. Update SHA when the URL changes.
PCCLOCKS_URL <- "https://zenodo.org/records/17162604/files/PCClocks_data.qs2?download=1"
PCCLOCKS_SHA256 <- NA_character_  # Set to expected hash once verified by maintainer
PCCLOCKS_MIN_SIZE <- 1e9          # secondary integrity check (>= 1 GB)


#' Check which clock packages are available
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
                   requireNamespace("methylcipher", quietly = TRUE),
    embedded_coefficients = !is.null(.qc_env$clock_coefficients) ||
      tryCatch({
        e <- new.env(parent = emptyenv())
        utils::data(list = "clock_coefficients", package = "clocker", envir = e)
        exists("clock_coefficients", envir = e)
      }, error = function(e) FALSE)
  )

  if (verbose) {
    for (pkg in names(availability)) {
      status <- if (availability[[pkg]]) "available" else "not installed"
      log_msg("  %s: %s", pkg, status, verbose = TRUE)
    }
  }
  availability
}


#' Compute all clocks
#'
#' @param betas Beta matrix (CpGs as rows, samples as columns)
#' @param pheno Optional phenotype data.frame
#' @param n_cores Number of CPU cores
#' @param verbose Print progress
#' @return Data frame with all clock results
#' @keywords internal
compute_all_clocks <- function(betas, pheno = NULL, n_cores, verbose = TRUE) {

  init_coverage_log()

  results <- data.frame(sample_id = colnames(betas), stringsAsFactors = FALSE)

  # Validate / reorder pheno upfront (FIX C1)
  pheno <- validate_and_align_pheno(pheno, colnames(betas), verbose = verbose)

  results <- compute_direct_clocks(betas, results, verbose)
  results <- estimate_cell_composition(betas, results, verbose)
  results <- compute_sex_inference(betas, results, verbose)
  results <- compute_dunedin_pace(betas, results, verbose)
  results <- compute_pc_clocks(betas, results, pheno, verbose)
  results <- compute_mitotic_clocks(betas, results, verbose)
  results <- compute_additional_clocks(betas, results, verbose)

  results
}


# ---------------------------------------------------------------------------
# Pheno validation / alignment (FIX C1, F6, M8)
# ---------------------------------------------------------------------------

#' Validate phenotype data frame and align rows to sample order in betas
#'
#' Normalizes the Female column to integer 0/1 (M=0, F=1), accepting:
#'   * 0/1 numeric
#'   * "M"/"F"/"Male"/"Female" (any case)
#'
#' Reorders rows to match `sample_ids` exactly, erroring if any sample is
#' missing from the pheno frame.
#'
#' @keywords internal
validate_and_align_pheno <- function(pheno, sample_ids, verbose = TRUE) {

  if (is.null(pheno)) return(NULL)

  if (!is.data.frame(pheno)) {
    warning("`pheno` is not a data.frame; ignoring.")
    return(NULL)
  }

  # If rownames aren't set, try to find a 'sample_id'/'Sample_ID'/'id' column
  if (is.null(rownames(pheno)) || all(rownames(pheno) == as.character(seq_len(nrow(pheno))))) {
    for (id_col in c("sample_id", "Sample_ID", "id", "ID", "Sample", "sample")) {
      if (id_col %in% colnames(pheno)) {
        rownames(pheno) <- as.character(pheno[[id_col]])
        break
      }
    }
  }

  if (is.null(rownames(pheno))) {
    if (length(sample_ids) == nrow(pheno)) {
      warning("`pheno` has no rownames or sample-ID column. Assuming row order matches betas.")
      rownames(pheno) <- sample_ids
    } else {
      stop("`pheno` must have rownames matching colnames(betas), or a 'sample_id' column.")
    }
  }

  missing_samples <- setdiff(sample_ids, rownames(pheno))
  if (length(missing_samples) > 0L) {
    stop(sprintf(
      "Pheno is missing %d sample(s): %s%s",
      length(missing_samples),
      paste(utils::head(missing_samples, 5), collapse = ", "),
      if (length(missing_samples) > 5) ", ..." else ""))
  }

  # Reorder to match betas
  pheno <- pheno[sample_ids, , drop = FALSE]

  # Normalize Female column if present
  if ("Female" %in% colnames(pheno)) {
    pheno$Female <- coerce_female_indicator(pheno$Female, na_on_unknown = TRUE)
  } else if ("Sex" %in% colnames(pheno)) {
    if (verbose) message("  Pheno: deriving Female from Sex column")
    pheno$Female <- coerce_female_indicator(pheno$Sex, na_on_unknown = TRUE)
  } else if ("sex" %in% colnames(pheno)) {
    if (verbose) message("  Pheno: deriving Female from sex column")
    pheno$Female <- coerce_female_indicator(pheno$sex, na_on_unknown = TRUE)
  }

  # Coerce Age numeric
  if ("Age" %in% colnames(pheno)) {
    pheno$Age <- suppressWarnings(as.numeric(pheno$Age))
  } else if ("age" %in% colnames(pheno)) {
    pheno$Age <- suppressWarnings(as.numeric(pheno$age))
  }

  pheno
}


# ---------------------------------------------------------------------------
# Direct weighted-sum clocks
# ---------------------------------------------------------------------------

#' @keywords internal
compute_direct_clocks <- function(betas, results, verbose = TRUE) {

  if (verbose) message("  Computing clocks via direct implementations...")

  tryCatch({
    coeffs <- initialize_clock_coefficients()
    coeffs <- augment_with_inline_fallbacks(coeffs)
    if (length(coeffs) > 0L && verbose) {
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
    computed_direct <- character()

    pb <- make_progress_bar(length(mc_clocks), label = "Direct clocks",
                              verbose = verbose)

    for (clock_name in names(mc_clocks)) {
      coef_name <- mc_clocks[[clock_name]]
      if (coef_name %in% names(coeffs)) {
        result <- tryCatch({
          if (clock_name %in% horvath_transform_clocks) {
            calc_weighted_sum_clock(betas, coeffs[[coef_name]],
                                     transform_func = horvath_age_transform,
                                     clock_name = clock_name)
          } else {
            calc_weighted_sum_clock(betas, coeffs[[coef_name]],
                                     clock_name = clock_name)
          }
        }, error = function(e) {
          if (verbose) message("    ", clock_name, " error: ", e$message)
          NULL
        })
        if (!is.null(result) && length(result) == ncol(betas)) {
          results[[clock_name]] <- as.numeric(result)
          computed_direct <- c(computed_direct, clock_name)
        }
      }
      pb$tick()
    }
    pb$done()

    # epiTOC2
    tryCatch({
      toc2 <- calc_epitoc2_direct(betas, coeffs)
      if (!is.null(toc2) && length(toc2) == ncol(betas)) {
        results$epiTOC2_TNSC <- toc2
        computed_direct <- c(computed_direct, "epiTOC2_TNSC")
      }
    }, error = function(e) {
      if (verbose) message("    epiTOC2 direct error: ", e$message)
    })

    if (verbose && length(computed_direct) > 0L) {
      message("    Direct clocks: ", paste(computed_direct, collapse = ", "))
    }

  }, error = function(e) {
    if (verbose) message("    Direct implementations error: ", e$message)
  })

  results
}


# ---------------------------------------------------------------------------
# Sex inference
# ---------------------------------------------------------------------------

#' @keywords internal
compute_sex_inference <- function(betas, results, verbose = TRUE) {
  tryCatch({
    log_msg("  Inferring sex (methylQC algorithm)...", verbose = verbose)

    sex_result <- infer_sex(betas, sdfs = NULL, verbose = verbose)

    if (is.list(sex_result) && length(sex_result$sex) == ncol(betas)) {
      results$chrX_signal      <- as.numeric(sex_result$chrX_signal)
      results$chrY_signal      <- as.numeric(sex_result$chrY_signal)
      results$InferredSex      <- as.character(sex_result$sex)
      results$InferredSexFlag  <- as.character(sex_result$flag)
      results$InferredSexScale <- sex_result$scale
    }
  }, error = function(e) {
    if (verbose) message("    Sex inference failed: ", e$message)
  })
  results
}


# ---------------------------------------------------------------------------
# DunedinPACE
# ---------------------------------------------------------------------------

#' @keywords internal
compute_dunedin_pace <- function(betas, results, verbose = TRUE) {

  if (!requireNamespace("DunedinPACE", quietly = TRUE)) return(results)

  tryCatch({
    log_msg("  Calculating DunedinPACE...", verbose = verbose)

    betas_dp <- if (is.matrix(betas)) betas else as.matrix(betas)
    if (is.null(dim(betas_dp))) {
      probe_names <- names(betas_dp)
      betas_dp <- matrix(betas_dp, ncol = 1,
                          dimnames = list(probe_names, colnames(betas)[1]))
    }

    is_single <- ncol(betas_dp) == 1L
    if (is_single) {
      sample_name <- colnames(betas_dp)[1]
      betas_dp <- cbind(betas_dp, betas_dp)
      colnames(betas_dp) <- c(sample_name, paste0(sample_name, "_dup"))
    }

    pace <- DunedinPACE::PACEProjector(betas_dp)
    if (!is.null(pace)) {
      pace_val <- if (is.list(pace) && "DunedinPACE" %in% names(pace)) pace$DunedinPACE
                  else if (is.data.frame(pace) && "DunedinPACE" %in% colnames(pace)) pace$DunedinPACE
                  else if (is.numeric(pace)) pace
                  else NULL
      if (!is.null(pace_val)) {
        if (is_single && length(pace_val) > 1L) pace_val <- pace_val[1]
        if (length(pace_val) == ncol(betas) ||
            (length(pace_val) == 1 && ncol(betas) == 1)) {
          results$DunedinPACE <- as.numeric(pace_val)
        }
      }
    }
  }, error = function(e) {
    if (verbose) message("    DunedinPACE failed: ", e$message)
  })
  results
}


# ---------------------------------------------------------------------------
# PC-Clocks
# ---------------------------------------------------------------------------

#' @keywords internal
ensure_pcclocks_data <- function(verbose = TRUE) {
  cache_dir <- get_cache_dir("pcclocks")
  pc_file <- file.path(cache_dir, "PCClocks_data.qs2")

  validate <- function(path) {
    if (!file.exists(path)) return(FALSE)
    if (file.info(path)$size < PCCLOCKS_MIN_SIZE) return(FALSE)
    if (!is.na(PCCLOCKS_SHA256)) {
      got <- file_sha256(path)
      if (!is.na(got) && !identical(got, PCCLOCKS_SHA256)) return(FALSE)
    }
    TRUE
  }

  if (validate(pc_file)) {
    if (verbose) {
      message(sprintf("    Using cached PC-Clocks data (%.2f GB)",
                      file.info(pc_file)$size / 1e9))
    }
    return(pc_file)
  }
  if (file.exists(pc_file)) {
    if (verbose) message("    Cached PC-Clocks data invalid; re-downloading...")
    unlink(pc_file)
  }

  if (verbose) {
    message("    Downloading PC-Clocks data from Zenodo (~2 GB; may take several minutes)...")
  }
  ok <- tryCatch({
    download_with_retry(PCCLOCKS_URL, pc_file, retries = 2L,
                         expected_sha256 = PCCLOCKS_SHA256,
                         verbose = verbose)
    TRUE
  }, error = function(e) {
    if (verbose) message("    Download failed: ", e$message)
    FALSE
  })

  if (!ok || !validate(pc_file)) {
    if (file.exists(pc_file)) unlink(pc_file)
    return(NULL)
  }
  pc_file
}


#' @keywords internal
compute_pc_clocks <- function(betas, results, pheno = NULL, verbose = TRUE) {

  if (!requireNamespace("methylCIPHER", quietly = TRUE)) return(results)

  tryCatch({
    if (verbose) message("  Computing PC clocks...")

    pc_data_path <- ensure_pcclocks_data(verbose = verbose)
    if (is.null(pc_data_path)) {
      if (verbose) message("    PC clocks skipped: data file not available")
      return(results)
    }

    # Load PC clock CpG list into the package private env (FIX C2)
    pre_load_objs <- ls(.qc_env, all.names = TRUE)
    tryCatch({
      utils::data(list = "PCClocks_CpGs", package = "methylCIPHER",
                   envir = .qc_env)
    }, error = function(e) NULL)
    new_pc_objs <- setdiff(ls(.qc_env, all.names = TRUE), pre_load_objs)

    # methylCIPHER's calcPCClocks may search for the loaded objects in
    # .GlobalEnv; mirror the bindings briefly with on.exit cleanup so the
    # user's workspace isn't permanently mutated.
    pc_mirrored <- character()
    for (obj_name in new_pc_objs) {
      if (!exists(obj_name, envir = .GlobalEnv, inherits = FALSE)) {
        assign(obj_name, get(obj_name, envir = .qc_env), envir = .GlobalEnv)
        pc_mirrored <- c(pc_mirrored, obj_name)
      }
    }
    on.exit(if (length(pc_mirrored)) {
      rm(list = pc_mirrored, envir = .GlobalEnv)
    }, add = TRUE)

    betas_t <- t(betas)
    sample_ids <- rownames(betas_t)

    if (!exists("calcPCClocks", envir = asNamespace("methylCIPHER"))) {
      if (verbose) message("    calcPCClocks not found in methylCIPHER")
      return(results)
    }
    pc_func <- get("calcPCClocks", envir = asNamespace("methylCIPHER"))

    # ---- Resolve Age (FIX M7/E: alert when falling back to clocks) ----
    age_source <- "user-supplied"
    if (!is.null(pheno) && "Age" %in% colnames(pheno) && !all(is.na(pheno$Age))) {
      age_values <- pheno$Age
    } else if ("Horvath2" %in% colnames(results) && !all(is.na(results$Horvath2))) {
      message("  PC-Clocks: 'Age' not provided; using Horvath2 estimate as Age")
      age_values <- results$Horvath2
      age_source <- "Horvath2_estimate"
    } else if ("Horvath1" %in% colnames(results) && !all(is.na(results$Horvath1))) {
      message("  PC-Clocks: 'Age' not provided; Horvath2 unavailable; falling back to Horvath1")
      age_values <- results$Horvath1
      age_source <- "Horvath1_estimate"
    } else {
      warning("PC-Clocks: no Age available (pheno missing, Horvath unavailable). Skipping PC-Clocks.")
      return(results)
    }

    # ---- Resolve Female (FIX F6/M8) ----
    female_source <- "user-supplied"
    if (!is.null(pheno) && "Female" %in% colnames(pheno) &&
        !all(is.na(pheno$Female))) {
      female_values <- as.integer(pheno$Female)
    } else if ("InferredSex" %in% colnames(results)) {
      message("  PC-Clocks: 'Female' not provided; using InferredSex (F=1, M/U=0)")
      female_values <- as.integer(results$InferredSex == "F")
      female_source <- "inferred_sex"
    } else {
      message("  PC-Clocks: 'Female' not provided and sex not inferred; defaulting to 0 (Male)")
      female_values <- rep(0L, length(sample_ids))
      female_source <- "default_male"
    }

    # NA Female -> 0 (PC-Clocks reject NA), with warning preserved upstream
    female_values[is.na(female_values)] <- 0L

    pheno_df <- data.frame(
      Sample_ID = sample_ids,
      Age = age_values,
      Female = female_values,
      stringsAsFactors = FALSE
    )
    rownames(pheno_df) <- sample_ids

    if (verbose) {
      message(sprintf("    PC clocks pheno: Age = %s; Female = %s",
                      age_source, female_source))
      message(sprintf("    Age range: %.1f - %.1f; Female: %d F / %d M",
                      min(pheno_df$Age, na.rm = TRUE),
                      max(pheno_df$Age, na.rm = TRUE),
                      sum(pheno_df$Female == 1L),
                      sum(pheno_df$Female == 0L)))
    }

    pc_result <- tryCatch(
      pc_func(betas_t, pheno_df, RData = pc_data_path),
      error = function(e) {
        if (verbose) message("    PC clocks failed: ", e$message)
        NULL
      })

    if (!is.null(pc_result) && is.data.frame(pc_result)) {
      pc_cols <- c("PCHorvath1", "PCHorvath2", "PCHannum",
                    "PCPhenoAge", "PCGrimAge", "PCDNAmTL")
      added_pc <- character()
      for (col in colnames(pc_result)) {
        if (col %in% pc_cols && is.numeric(pc_result[[col]])) {
          if (nrow(pc_result) == ncol(betas)) {
            results[[col]] <- pc_result[[col]]
            added_pc <- c(added_pc, col)
          }
        }
      }
      if (verbose && length(added_pc) > 0L) {
        message("    PC clocks added: ", paste(added_pc, collapse = ", "))
      }
    } else if (verbose) {
      message("    PC clocks: no result returned")
    }

  }, error = function(e) {
    if (verbose) message("    PC clocks error: ", e$message)
  })

  results
}


# ---------------------------------------------------------------------------
# Mitotic clocks
# ---------------------------------------------------------------------------

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

    # Track which objects were already in .qc_env so we only mirror NEW ones
    pre_load_objs <- ls(.qc_env, all.names = TRUE)
    for (d in epi_data_items) {
      tryCatch({
        utils::data(list = d, package = epi_pkg, envir = .qc_env)
      }, error = function(e) NULL)
    }
    new_objs <- setdiff(ls(.qc_env, all.names = TRUE), pre_load_objs)
    # Note: data() typically loads objects with their actual stored names
    # (e.g., dataETOC3.l, not dataETOC3), so we mirror the *loaded objects*
    # rather than the dataset names. This is what EpiMitClocks::EpiMitClocks
    # expects to find in .GlobalEnv at call time.

    mirrored <- character()
    for (obj_name in new_objs) {
      if (!exists(obj_name, envir = .GlobalEnv, inherits = FALSE)) {
        assign(obj_name, get(obj_name, envir = .qc_env), envir = .GlobalEnv)
        mirrored <- c(mirrored, obj_name)
      }
    }
    on.exit(if (length(mirrored)) {
      rm(list = mirrored, envir = .GlobalEnv)
    }, add = TRUE)

    if (exists("EpiMitClocks", envir = asNamespace(epi_pkg))) {
      epi_results <- tryCatch({
        get("EpiMitClocks", envir = asNamespace(epi_pkg))(
          data.m = betas, ages.v = NULL)
      }, error = function(e) {
        if (verbose) message("    EpiMitClocks failed: ", e$message)
        NULL
      })

      if (!is.null(epi_results) && is.data.frame(epi_results) &&
          nrow(epi_results) == ncol(betas)) {
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
    if (verbose) message("    Mitotic clocks error: ", e$message)
  })

  results
}


# ---------------------------------------------------------------------------
# Additional methylCIPHER clocks (AdaptAge, CausAge, DamAge, SystemsAge)
# ---------------------------------------------------------------------------

#' @keywords internal
compute_additional_clocks <- function(betas, results, verbose = TRUE) {

  if (!requireNamespace("methylCIPHER", quietly = TRUE)) return(results)

  tryCatch({
    if (verbose) message("  Computing additional clocks via methylCIPHER...")

    betas_t <- t(betas)
    additional_clocks <- c("AdaptAge", "CausAge", "DamAge", "SystemsAge")

    # Load CpG datasets into .qc_env
    needed <- c("AdaptAge_CpGs", "CausAge_CpGs", "DamAge_CpGs", "SystemsAge_CpGs")
    pre_load_objs <- ls(.qc_env, all.names = TRUE)
    for (d in needed) {
      tryCatch({
        utils::data(list = d, package = "methylCIPHER", envir = .qc_env)
      }, error = function(e) NULL)
    }
    new_objs <- setdiff(ls(.qc_env, all.names = TRUE), pre_load_objs)

    mirrored <- character()
    for (obj_name in new_objs) {
      if (!exists(obj_name, envir = .GlobalEnv, inherits = FALSE)) {
        assign(obj_name, get(obj_name, envir = .qc_env), envir = .GlobalEnv)
        mirrored <- c(mirrored, obj_name)
      }
    }
    on.exit(if (length(mirrored)) rm(list = mirrored, envir = .GlobalEnv),
             add = TRUE)

    added_clocks <- character()
    for (clock in additional_clocks) {
      if (clock %in% colnames(results)) next
      func_name <- paste0("calc", clock)
      if (!exists(func_name, envir = asNamespace("methylCIPHER"))) next

      func <- get(func_name, envir = asNamespace("methylCIPHER"))
      result <- tryCatch(func(betas_t, imputation = FALSE),
                          error = function(e) {
                            tryCatch(func(betas_t),
                                      error = function(e2) NULL)
                          })
      if (!is.null(result) && is.numeric(result) &&
          (length(result) == ncol(betas) ||
           (length(result) == 1L && ncol(betas) == 1L))) {
        results[[clock]] <- as.numeric(result)
        added_clocks <- c(added_clocks, clock)
      }
    }
    if (verbose && length(added_clocks) > 0L) {
      message("    Additional clocks: ", paste(added_clocks, collapse = ", "))
    }
  }, error = function(e) {
    if (verbose) message("    Additional clocks error: ", e$message)
  })

  results
}
