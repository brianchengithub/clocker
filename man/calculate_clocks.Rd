\name{calculate_clocks}
\alias{calculate_clocks}
\title{Calculate Epigenetic Clocks}
\description{
  Unified function for calculating 40+ epigenetic clocks from DNA methylation data.
  Automatically detects input type (IDAT directory or beta matrix) and array platform.
  Includes EpiDISH cell type deconvolution, PC-adjusted clocks, and age acceleration metrics.
}
\usage{
calculate_clocks(input, n_cores = NULL, verbose = TRUE)
}
\arguments{
  \item{input}{Either a character string path to an IDAT directory, or a numeric matrix 
    of beta values with CpG probe IDs as row names and sample IDs as column names.}
  \item{n_cores}{Number of CPU cores for parallel processing. NULL (default) auto-detects 
    the optimal number based on data size and available resources.}
  \item{verbose}{Logical. Print progress messages during calculation. Default TRUE.}
}
\value{
  A \code{data.frame} with one row per sample and columns including:
  \itemize{
    \item \code{Sample_ID} - Sample identifiers
    \item \code{CellType_RPC_*} - EpiDISH RPC cell type estimates (7 types)
    \item \code{CellType_CP_*} - EpiDISH CP cell type estimates (7 types)
    \item \code{InferredSex} - Chromosome-based sex inference (M/F/U)
    \item \code{Horvath1}, \code{Horvath2}, \code{Hannum}, \code{PhenoAge}, etc. - Age clocks
    \item \code{PCHorvath1}, \code{PCHorvath2}, etc. - PC-adjusted clocks
    \item \code{DunedinPACE} - Pace of biological aging
    \item \code{epiTOC2_TNSC} - Mitotic clock
  }
}
\details{
  When provided IDAT files, preprocessing uses SeSAMe's openSesame() pipeline 
  (NOOB normalization, dye-bias correction, pOOBAH).

  Cell type deconvolution uses EpiDISH with the centDHSbloodDMC.m reference matrix, 
  providing estimates for 7 blood cell types (B, CD4T, CD8T, Mono, Neu, NK, Eos) 
  via both RPC (Robust Partial Correlations) and CP (Constrained Projection) methods.

  PC clocks require a ~2GB data file automatically downloaded from Zenodo on first use.
  When Age and Sex phenotype data are not provided, PC clocks use Horvath2 age estimates 
  and chromosome-based sex inference.

  Missing CpG values are imputed using reference median values from the included 
  reference dataset.
}
\examples{
\dontrun{
  # From IDAT directory
  results <- calculate_clocks("/path/to/idats/")

  # From beta matrix
  results <- calculate_clocks(beta_matrix)

  # With manual core count
  results <- calculate_clocks(beta_matrix, n_cores = 4)
}
}
\seealso{
  \code{\link[EpiDISH]{epidish}}, \code{\link[sesame]{openSesame}}
}
