Package: quickclocks
Type: Package
Title: Unified Epigenetic Clock Calculator
Version: 2.0.0
Date: 2025-01-09
Authors@R: c(
    person("Brian", "Chen", email = "brian.chen4@ucsf.edu", role = c("aut", "cre"))
  )
Description: A simple, unified interface for calculating 40+ epigenetic clocks from 
    DNA methylation data. Automatically detects input type (IDAT directory or beta 
    matrix) and array platform (MSA, EPIC+, EPIC, 450K, 27K). Integrates clocks from 
    SeSAMe, DunedinPACE, PC-Clocks, EpiMitClocks (epiTOC2), and methylCIPHER packages.
    Includes EpiDISH cell type deconvolution (RPC and CP methods) and smart thread 
    optimization based on data size and available RAM.
    Returns results as a clean data.frame with sample IDs and all clock values.
License: MIT + file LICENSE
URL: https://github.com/brianchengithub/quickclocks
BugReports: https://github.com/brianchengithub/quickclocks/issues
Encoding: UTF-8
LazyData: true
RoxygenNote: 7.3.1
Depends:
    R (>= 4.0.0)
Imports:
    parallel,
    doParallel,
    foreach,
    data.table,
    matrixStats,
    methods,
    stats,
    utils
Suggests:
    sesame,
    sesameData,
    EpiDISH,
    DunedinPACE,
    PCClocks,
    EpiMitClocks,
    methylCIPHER,
    qs2,
    minfi,
    IlluminaHumanMethylation27kmanifest,
    IlluminaHumanMethylation450kmanifest,
    IlluminaHumanMethylationEPICmanifest,
    IlluminaHumanMethylationEPICv2manifest,
    impute,
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
biocViews: 
    DNAMethylation,
    Epigenetics,
    Preprocessing,
    QualityControl,
    MethylationArray
VignetteBuilder: knitr
Config/testthat/edition: 3
