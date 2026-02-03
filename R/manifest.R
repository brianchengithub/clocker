# ============================================================================
# Manifest Download and Caching Functions
# Downloads and caches Illumina array manifests for sex inference
# ============================================================================


#' Get manifest cache directory
#' @return Path to cache directory
#' @keywords internal
get_manifest_cache_dir <- function() {
  cache_dir <- file.path(Sys.getenv("HOME"), ".epigenetic_clock_calculator", "manifests")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  return(cache_dir)
}


#' Get manifest URL for a platform
#' @param platform Platform name (EPIC, EPICv2, HM450, etc.)
#' @return URL to manifest file, or NULL if unknown
#' @keywords internal
get_manifest_url <- function(platform) {
  base_url <- "https://github.com/zhou-lab/InfiniumAnnotationV1/raw/main/Anno"
  urls <- list(
    "MSA"   = paste0(base_url, "/MSA/MSA.hg38.manifest.tsv.gz"),
    "EPICv2" = paste0(base_url, "/EPICv2/EPICv2.hg38.manifest.tsv.gz"),
    "EPIC+"  = paste0(base_url, "/EPIC+/EPIC+.hg38.manifest.tsv.gz"),
    "EPIC"   = paste0(base_url, "/EPIC/EPIC.hg38.manifest.tsv.gz"),
    "HM450"  = paste0(base_url, "/HM450/HM450.hg38.manifest.tsv.gz")
  )
  return(urls[[platform]])
}


#' Download and cache manifest for a platform
#' @param platform Platform name
#' @param verbose Print progress
#' @return Data frame with manifest data, or NULL on failure
#' @keywords internal
download_manifest <- function(platform, verbose = TRUE) {
  cache_dir <- get_manifest_cache_dir()

  use_qs2 <- requireNamespace("qs2", quietly = TRUE)
  cache_file_qs2 <- file.path(cache_dir, paste0(platform, ".hg38.manifest.qs2"))
  cache_file_rds <- file.path(cache_dir, paste0(platform, ".hg38.manifest.rds"))

  # Check for cached file (prefer qs2, then rds)
  if (use_qs2 && file.exists(cache_file_qs2)) {
    if (verbose) message("    Loading cached manifest for ", platform, " (qs2)")
    return(qs2::qs_read(cache_file_qs2))
  } else if (file.exists(cache_file_rds)) {
    if (verbose) message("    Loading cached manifest for ", platform, " (rds)")
    manifest <- readRDS(cache_file_rds)
    if (use_qs2) {
      tryCatch({
        qs2::qs_save(manifest, cache_file_qs2)
        unlink(cache_file_rds)
        if (verbose) message("    Upgraded cache to qs2 format")
      }, error = function(e) NULL)
    }
    return(manifest)
  }

  # Download
  url <- get_manifest_url(platform)
  if (is.null(url)) {
    if (verbose) message("    Unknown platform: ", platform)
    return(NULL)
  }

  if (verbose) message("    Downloading manifest for ", platform, "...")

  temp_file <- tempfile(fileext = ".tsv.gz")

  tryCatch({
    download.file(url, temp_file, mode = "wb", quiet = !verbose)

    manifest <- read.table(gzfile(temp_file), header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, quote = "",
                          comment.char = "", fill = TRUE)

    if (use_qs2) {
      qs2::qs_save(manifest, cache_file_qs2)
      if (verbose) message("    Cached manifest: ", nrow(manifest), " probes (qs2 format)")
    } else {
      saveRDS(manifest, cache_file_rds)
      if (verbose) message("    Cached manifest: ", nrow(manifest), " probes (rds format)")
    }

    unlink(temp_file)
    return(manifest)

  }, error = function(e) {
    if (verbose) message("    Failed to download manifest: ", e$message)
    unlink(temp_file)
    return(NULL)
  })
}


#' Get chromosome information for probes from manifest
#' @param platform Platform name
#' @param verbose Print progress
#' @return Named vector: probe_id -> chromosome
#' @keywords internal
get_probe_chromosomes <- function(platform, verbose = TRUE) {
  manifest <- download_manifest(platform, verbose)

  if (is.null(manifest)) return(NULL)

  # Find probe ID column
  probe_col <- NULL
  for (col in c("Probe_ID", "probeID", "IlmnID", "Name")) {
    if (col %in% colnames(manifest)) {
      probe_col <- col
      break
    }
  }
  if (is.null(probe_col)) probe_col <- colnames(manifest)[1]

  # Find chromosome column
  chr_col <- NULL
  for (col in c("CpG_chrm", "CHR", "chr", "Chromosome", "seqnames")) {
    if (col %in% colnames(manifest)) {
      chr_col <- col
      break
    }
  }

  if (is.null(chr_col)) {
    if (verbose) message("    Could not find chromosome column in manifest")
    return(NULL)
  }

  probe_chr <- manifest[[chr_col]]
  names(probe_chr) <- manifest[[probe_col]]

  return(probe_chr)
}
