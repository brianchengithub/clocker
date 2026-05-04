# ============================================================================
# Cell Type Deconvolution
# EpiDISH-based cell composition (RPC and CP methods).
# ============================================================================


#' Estimate cell type composition using EpiDISH
#'
#' Runs RPC (robust partial correlations) and CP (constrained projection)
#' against the centDHSbloodDMC.m reference (7 immune cell types).
#'
#' @keywords internal
estimate_cell_composition <- function(betas, results, verbose = TRUE) {

  if (!requireNamespace("EpiDISH", quietly = TRUE)) {
    if (verbose) {
      message("  EpiDISH not installed - skipping cell deconvolution")
      message("  Install with: BiocManager::install('EpiDISH')")
    }
    return(results)
  }

  tryCatch({
    log_msg("  Estimating cell composition with EpiDISH...", verbose = verbose)

    ref_matrix <- NULL
    ref_name <- NULL
    for (ref in c("centDHSbloodDMC.m")) {
      tryCatch({
        e <- new.env(parent = emptyenv())
        utils::data(list = ref, package = "EpiDISH", envir = e)
        if (exists(ref, envir = e)) {
          ref_matrix <- get(ref, envir = e)
          ref_name <- ref
          break
        }
      }, error = function(e) NULL)
    }
    if (is.null(ref_matrix)) {
      if (verbose) message("    Could not load EpiDISH reference matrix")
      return(results)
    }

    if (verbose) {
      message("    Loaded reference: ", ref_name,
              " (", ncol(ref_matrix), " cell types)")
      message("    Cell types: ",
              paste(colnames(ref_matrix), collapse = ", "))
    }

    common_probes <- intersect(rownames(betas), rownames(ref_matrix))
    if (verbose) {
      message("    Overlapping probes with reference: ", length(common_probes))
    }
    if (length(common_probes) < 100) {
      if (verbose) {
        message("    Warning: very few overlapping probes; ",
                "estimates may be unreliable")
      }
    }
    if (length(common_probes) == 0L) {
      if (verbose) message("    No overlap with reference; skipping")
      return(results)
    }

    betas_subset <- betas[common_probes, , drop = FALSE]
    ref_subset   <- ref_matrix[common_probes, , drop = FALSE]

    # RPC method
    rpc_result <- tryCatch(
      EpiDISH::epidish(beta.m = betas_subset, ref.m = ref_subset, method = "RPC"),
      error = function(e) {
        if (verbose) message("    EpiDISH RPC error: ", e$message)
        NULL
      })
    if (!is.null(rpc_result) && !is.null(rpc_result$estF)) {
      rpc_fracs <- as.data.frame(rpc_result$estF)
      if (nrow(rpc_fracs) == ncol(betas)) {
        colnames(rpc_fracs) <- paste0("CellType_RPC_", colnames(rpc_fracs))
        for (col in colnames(rpc_fracs)) results[[col]] <- rpc_fracs[[col]]
        if (verbose) {
          message("    EpiDISH RPC: ", ncol(rpc_fracs), " cell types estimated")
        }
      }
    }

    # CP method
    cp_result <- tryCatch(
      EpiDISH::epidish(beta.m = betas_subset, ref.m = ref_subset, method = "CP"),
      error = function(e) {
        if (verbose) message("    EpiDISH CP error: ", e$message)
        NULL
      })
    if (!is.null(cp_result) && !is.null(cp_result$estF)) {
      cp_fracs <- as.data.frame(cp_result$estF)
      if (nrow(cp_fracs) == ncol(betas)) {
        colnames(cp_fracs) <- paste0("CellType_CP_", colnames(cp_fracs))
        for (col in colnames(cp_fracs)) results[[col]] <- cp_fracs[[col]]
        if (verbose) {
          message("    EpiDISH CP: ", ncol(cp_fracs), " cell types estimated")
        }
      }
    }

  }, error = function(e) {
    if (verbose) message("    EpiDISH deconvolution error: ", e$message)
  })

  results
}
