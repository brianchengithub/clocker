\name{quickclocks-package}
\alias{quickclocks-package}
\alias{quickclocks}
\docType{package}
\title{Unified Epigenetic Clock Calculator}
\description{
  A comprehensive R package for calculating 40+ epigenetic clocks from DNA methylation 
  data. Integrates clocks from multiple research labs into a single, easy-to-use 
  interface with a single function call.

  Key features:
  \itemize{
    \item 40+ epigenetic clocks from multiple research groups
    \item Single function interface - one call calculates everything
    \item Automatic input detection (IDAT files or beta matrices)
    \item EpiDISH cell type deconvolution (RPC and CP methods)
    \item PC-adjusted clocks with automatic phenotype inference
    \item Parallel processing with automatic core detection
    \item Smart caching for manifests and large data files
  }
}
\details{
  The main function is \code{\link{calculate_clocks}}, which accepts either 
  a directory of IDAT files or a beta value matrix and returns a data.frame 
  with all calculated clocks.

  Supported platforms: EPICv2, EPIC, 450K, 27K, MSA.
}
\author{
  Author Name
}
\seealso{
  \code{\link{calculate_clocks}}
}
