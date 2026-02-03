\name{calculate_clocks}
\alias{calculate_clocks}
\alias{run_epigenetic_clocks}
\title{Calculate Epigenetic Clocks}
\description{
Main function to calculate 40+ epigenetic clocks from DNA methylation data.
Automatically detects whether input is a directory path (IDAT files) or
a numeric matrix of beta values.
}
\usage{
calculate_clocks(input, pheno = NULL, n_cores = NULL, verbose = TRUE)

run_epigenetic_clocks(input, pheno = NULL, n_cores = NULL, verbose = TRUE)
}
\arguments{
  \item{input}{Either a character string path to a directory containing IDAT files,
    a path to an RDS or CSV file, or a numeric matrix of beta values with CpG
    probe IDs as rownames and sample IDs as colnames.}
  \item{pheno}{Optional data.frame with phenotype data. Should have row names
    matching sample IDs, an \code{Age} column (numeric), and a \code{Female}
    column (1=female, 0=male, 0.5=unknown). Required for accurate PC clock
    calculations.}
  \item{n_cores}{Number of CPU cores to use. Default \code{NULL} selects
    automatically based on data size, available cores, and RAM.}
  \item{verbose}{Logical. Print progress messages. Default \code{TRUE}.}
}
\value{
A data.frame with sample IDs as rows and columns including:
  \itemize{
    \item \code{sample_id}: Sample identifier
    \item \code{Horvath1}, \code{Horvath2}, \code{Hannum}, etc.: Clock values
    \item \code{CellType_RPC_*}: EpiDISH RPC cell fractions
    \item \code{CellType_CP_*}: EpiDISH CP cell fractions
    \item \code{InferredSex}: Predicted sex (F/M/U)
    \item \code{DunedinPACE}: Pace of aging
    \item \code{PCHorvath1}, \code{PCPhenoAge}, etc.: PC clock values
  }
}
\details{
\strong{Preprocessing:} IDAT files are preprocessed with SeSAMe using NOOB
normalization, dye-bias correction, and pOOBAH masking.

\strong{Cell Deconvolution:} EpiDISH estimates 7 blood cell types (B, CD4T,
CD8T, Mono, Neu, NK, Eos) using both RPC and CP methods.

\strong{PC Clocks:} Require ~2GB data file from Zenodo (auto-downloaded and
cached on first use). Age and sex are auto-inferred if not provided.

\strong{Imputation:} Missing values are imputed using reference betas
(if available) or row medians.

\code{run_epigenetic_clocks} is an alias for backward compatibility.
}
\examples{
\dontrun{
# From IDAT directory
results <- calculate_clocks("/path/to/idat/directory")

# From beta matrix
results <- calculate_clocks(my_beta_matrix)

# With phenotype data for PC clocks
pheno <- data.frame(Age = c(45, 52), Female = c(1, 0),
                    row.names = colnames(my_betas))
results <- calculate_clocks(my_betas, pheno = pheno)
}
}
\seealso{
\code{\link[EpiDISH]{epidish}}, \code{\link[sesame]{openSesame}}
}
