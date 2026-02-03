\name{quickclocks-package}
\alias{quickclocks-package}
\alias{quickclocks}
\docType{package}
\title{Unified Epigenetic Clock Calculator}
\description{
Calculate 40+ epigenetic clocks from DNA methylation data with a single
function call. Integrates SeSAMe for IDAT preprocessing, EpiDISH for cell
type deconvolution, DunedinPACE, PC-Clocks, epiTOC2, and methylCIPHER.

Supports EPICv2/EPIC+, EPIC, 450K, 27K, and MSA arrays.
}
\details{
The main function is \code{\link{calculate_clocks}}, which accepts either
a directory of IDAT files or a beta value matrix and returns a data.frame
with all available clock results.
}
\author{Brian Chen}
\seealso{
\code{\link{calculate_clocks}} for the main entry point.
}
