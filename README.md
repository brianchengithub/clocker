# clocker

**Unified Epigenetic Clock Calculator**

A comprehensive R package for calculating 40+ epigenetic clocks from DNA
methylation data with a single function call. Integrates clocks from
multiple research labs into one interface, accepts IDAT files or beta
matrices, and handles imputation, sex inference, cell deconvolution, and
clock projection automatically.

## Features

- **40+ epigenetic clocks** from across the literature in one place
- **Single function call**: `clocker(input)` does everything
- **Many input formats**: IDAT directory, IDAT file paths, beta matrix,
  or `.rds` / `.qs2` / `.csv` / `.tsv` file
- **EPICv2 native**: probe-name normalization with EPIC/450K position
  matching where possible
- **kNN imputation** with reference-mean zero-shot fallback
- **Sex inference** (exact algorithm, verified bit-for-bit)
- **EpiDISH cell deconvolution** (RPC + CP, 7 blood cell types)
- **PC-Clocks** with automatic Age/Female fallback when pheno is missing
- **Embedded coefficients** for the small clocks — no external R-package
  dependency required for Horvath, Hannum, PhenoAge, DNAmTL, Lin, Zhang,
  Bocklandt, Weidner, Vidal-Bralo, Garagnani, AdaptAge, CausAge, DamAge,
  SystemsAge, epiTOC2 once the build script has run
- **CRAN-compliant cache**: `tools::R_user_dir("clocker", "cache")`
- **Per-sample missing-probe report** to spot samples where any clock's
  coverage is degraded

## Supported Platforms

| Platform        | Probes      | Status            |
| --------------- | ----------- | ----------------- |
| EPICv2 / EPIC+  | ~930,000    | Fully supported   |
| EPIC v1         | ~865,000    | Fully supported   |
| 450K            | ~485,000    | Fully supported   |
| 27K             | ~27,000     | Supported         |
| MSA             | ~285,000    | Supported         |

## Installation

Install dependencies first (one-time setup):

```r
source("install_dependencies.R")
```

Then install clocker from GitHub:

```r
install.packages("pak")           # one-time
pak::pak("brianchengithub/clocker")
```

The installer pulls the required Bioconductor (`sesame`, `EpiDISH`) and
GitHub-only (`DunedinPACE`, `methylCIPHER`, `EpiMitClocks`) packages.
Optional helpers (`qs2`, `digest`, `progress`, `matrixStats`, `curl`) are
also installed; clocker has graceful fallbacks if any are missing.

## Quick start

```r
library(clocker)

# All clocks from an IDAT directory
results <- clocker("/filepath/idats/")

# Or from a pre-processed beta matrix (CpGs x samples, [0, 1] values)
results <- clocker(beta_matrix)

# With phenotype info — pheno$Age and pheno$Sex auto-coerced
pheno <- data.frame(
  Age = c(45, 67, 22),
  Sex = c("Female", "M", "f"),       # accepts 0/1, M/F, Male/Female (any case)
  row.names = c("S1", "S2", "S3")
)
results <- clocker(beta_matrix, pheno = pheno)
```

### Common options

```r
results <- clocker(
  input,
  pheno                   = NULL,    # optional phenotype data.frame
  n_cores                 = NULL,    # NULL = auto-detect
  reference_path          = NULL,    # path to reference_betas.rds
  knn_k                   = 10,      # neighbours for kNN imputation
  knn_zero_shot_threshold = 0.10,    # > this fraction missing -> zero-shot
  missing_report_path     = NULL,    # CSV path for per-sample coverage
  verbose                 = TRUE
)
```

## Output

A `data.frame` with one row per sample. Columns are grouped:

| Prefix / Column          | Description                                          |
| ------------------------ | ---------------------------------------------------- |
| `sample_id`              | Sample identifier                                    |
| `InferredSex`            | methylQC call: `F`, `M`, or `Unclear`                |
| `InferredSexFlag`        | `F_band_only` / `M_band_only` / `outside_both_bands` / etc. |
| `InferredSexScale`       | `intensity` (IDAT input) or `beta` (matrix input)    |
| `chrX_signal`, `chrY_signal` | Per-sample chrX/chrY summary stat used by sex caller |
| `Clock_*`                | First-generation clocks (Horvath1, Horvath2, Hannum, PhenoAge, DNAmTL, Lin, Zhang, ...) |
| `PC_*`                   | PC-adjusted clocks (PCHorvath1/2, PCHannum, PCPhenoAge, PCGrimAge, PCDNAmTL) |
| `DunedinPACE`            | Pace-of-aging estimate                               |
| `epiTOC2_TNSC`           | Total stem cell divisions                            |
| `EpiCMIT`, `MiAge`, `Replitali` | Other mitotic clocks                          |
| `AdaptAge`, `CausAge`, `DamAge`, `SystemsAge` | causality / damage clocks      |
| `CellType_RPC_*`         | EpiDISH RPC fractions (7 blood cell types)           |
| `CellType_CP_*`          | EpiDISH CP fractions                                 |
| `Accel_*`                | Age acceleration residuals (when `pheno$Age` provided) |

Per-sample imputation diagnostics are stashed in
`clocker:::.qc_env$imputation_info` after every run.

## How it works (high level)

1. **Load + validate.** IDAT files go through SeSAMe (`openSesame(prep="QCDPB")`).
   Beta matrices are sanity-checked and clipped to `[0, 1]`. Robust
   M-value detection (median + IQR) auto-converts M-values when needed.
2. **EPICv2 harmonization.** Suffix-bearing probe names (`cg00000029_TC11`)
   are mapped to base CpG IDs. When multiple replicates exist, the one
   matching the original EPIC v1 / 450K probe by genomic position +
   design + strand is preferred over averaging.
3. **kNN imputation** with reference-mean zero-shot fallback. See [Imputation](#imputation) below.
4. **Sex inference.** Exact port of methylQC's algorithm —
   curated chrX/chrY probe panel, threshold sweep, cluster regression,
   ±5σ orthogonal-distance bands. Bit-for-bit verified against the
   methylQC reference. See [Sex inference](#sex-inference) below.
5. **Clock projection.** All available clocks are computed; PC-Clocks
   data (~2 GB) downloads on first use to the cache.
6. **Cell deconvolution.** EpiDISH RPC + CP against `centDHSbloodDMC.m`.
7. **Output assembly.** Stable column ordering with `Clock_*`, `PC_*`,
   `Accel_*`, `CellType_*` prefixes.

## Imputation

Two-pass hybrid:

- **Within-batch kNN.** For each sample with < `knn_zero_shot_threshold`
  missing-probe fraction, the `k` most similar OTHER samples in the input
  batch are found by `1 - Pearson` distance over probes both samples have
  non-missing. Each missing probe is imputed as a Gaussian-kernel-weighted
  average across those neighbours.
- **Zero-shot reference fallback.** Samples above the threshold (or any
  residual probes the kNN neighbours also lack, or single-sample runs)
  fall through to the per-probe means in your `reference_betas.rds`.

Reference resolution order:

1. `reference_path` argument to `clocker()`
2. `options(clocker.reference_betas = "/path/...")`
3. `Sys.getenv("CLOCKER_REFERENCE_BETAS")`
4. `system.file("extdata", "reference_betas.rds", package = "clocker")`

`reference_betas.rds` should be a **named numeric vector** of mean beta
values keyed by CpG ID. Data-frame and matrix forms are also accepted
for backward compat.

## Sex inference

Exact port of the methylQC algorithm:

1. **Per-sample sex signal.** When SigDFs are available (IDAT input),
   total signal intensity (`MG + MR + UG + UR`) is computed over a
   curated panel of 314 chrY probes (PAR-excluded, cross-hyb removed)
   and 3,433 chrX-X-inactivation probes (PAR-excluded). For beta-only
   input, median beta over the same probes substitutes.
2. **Threshold optimization.** Sweep candidate chrY thresholds in
   `quantile(chrY, seq(0.15, 0.85, 0.01))` and pick the one minimizing
   total absolute residuals from per-cluster `lm(chrY ~ chrX)` fits.
3. **Confidence bands.** Refit each cluster, compute orthogonal distance
   from each sample to each line, set band thresholds at 5σ.
4. **Classification.** F-band only → F; M-band only → M; both bands →
   tiebreak by chrY threshold; neither band → Unclear.

For batches < 6 samples (data-driven fit impossible), embedded reference
parameters in `DEFAULT_SEX_REFERENCE` are used. Recalibrate for atypical
tissues (tumor, placenta, brain) via:

```r
set_sex_reference(
  scale     = "intensity",                # or "beta"
  threshold = 6500,
  male      = list(slope = -0.10, intercept = 7800, sigma = 350),
  female    = list(slope = -0.02, intercept = 1100, sigma = 300)
)
```

## Caching

clocker caches large files in CRAN-compliant locations resolved by:

1. `getOption("clocker.cache_dir")` if set
2. `Sys.getenv("CLOCKER_CACHE")` if set
3. `tools::R_user_dir("clocker", "cache")` (XDG-compliant on Linux)
4. `tempdir()/clocker` as a final fallback

Cached items:

- `manifests/EPIC.hg38.manifest.qs2` (and platform variants)
- `pcclocks/PCClocks_data.qs2` (~2 GB; SHA-256 verified when set)

## Performance

- ~20 seconds for 60+ EPIC samples after first run (manifests cached)
- First run downloads array manifests (~70 MB) and PC-Clocks data (~2 GB)
- Memory ≈ 5 GB per 100 EPICv2 samples (heavy phase is PC-Clocks projection)
- `n_cores` auto-detection accounts for available RAM and dataset size

## Troubleshooting

**`EpiDISH not installed`** — `BiocManager::install("EpiDISH")`

**PC-Clocks data download hangs or times out** — Manually download
`PCClocks_data.qs2` from
[Zenodo](https://zenodo.org/records/13952402) and place it at
`tools::R_user_dir("clocker", "cache")/pcclocks/PCClocks_data.qs2`.

**Low probe-overlap warnings** — Check that input rownames are standard
Illumina IDs (e.g., `cg00000029`). On EPICv2 data the package strips
suffixes automatically; if you see this warning on EPICv2 input, the
underlying data may have been renamed before reaching clocker.

**Diagnostic dump** — `Rscript diagnose_packages.R > diag.log 2>&1` and
share `diag.log` when reporting issues.

## File structure

```
clocker/
├── DESCRIPTION
├── NAMESPACE
├── README.md
├── METHODS.md
├── install_dependencies.R       # One-time dependency setup
├── diagnose_packages.R          # Troubleshooting utility
├── R/
│   ├── clocker.R                # Main function: clocker() + aliases
│   ├── input_processing.R       # IDAT loading, betas validation
│   ├── imputation.R             # kNN + zero-shot vs. reference_betas.rds
│   ├── sex_inference.R          # methylQC algorithm (verbatim)
│   ├── sex_probe_lists.R        # 314 chrY + 3,433 chrX probes (MIT, attribution)
│   ├── manifest.R               # Array manifest download/cache
│   ├── clock_engines.R          # Weighted-sum + epiTOC2 + missing-probe report
│   ├── clock_coefficients.R     # Embedded coefficients + live fallbacks
│   ├── clock_computation.R      # Orchestrator: routes to each clock backend
│   ├── cell_deconvolution.R     # EpiDISH RPC + CP
│   ├── utils.R                  # Logging, validation, cache, F/M coercion
│   └── zzz.R                    # Package load hooks + private env
├── data-raw/
│   └── build_coefficients.R     # Maintainer-only: builds data/clock_coefficients.rda
├── data/
│   └── clock_coefficients.rda   # Generated (run data-raw script once)
└── inst/
    └── extdata/
        └── reference_betas.rds  # Probe-mean reference for zero-shot imputation
```

## Citation

If you use clocker, please cite the original publications:

**Core methods**

- **SeSAMe**: Zhou W et al. (2018) *Nucleic Acids Research* 46:e123
- **EpiDISH**: Teschendorff AE et al. (2017) *BMC Bioinformatics* 18:105
- **DunedinPACE**: Belsky DW et al. (2022) *eLife* 11:e73420
- **PC-Clocks**: Higgins-Chen AT et al. (2022) *Nature Aging* 2:644
- **epiTOC2**: Teschendorff AE (2020) *Genome Medicine* 12:56
- **methylCIPHER**: Higgins-Chen AT et al. (2022)
- **methylQC** (sex caller): Cheng B (https://github.com/brianchengithub/methylQC), MIT

**Original clocks**

- **Horvath**: Horvath S (2013) *Genome Biology* 14:R115; Horvath et al. (2018) *Aging* 10:1758
- **Hannum**: Hannum G et al. (2013) *Molecular Cell* 49:359
- **PhenoAge**: Levine ME et al. (2018) *Aging* 10:573
- **GrimAge**: Lu AT et al. (2019) *Aging* 11:303
- **DNAmTL**: Lu AT et al. (2019) *Aging* 11:5895

## License

MIT — see LICENSE.

## Contributing

Issues and pull requests welcome.

## Changelog

See `CHANGES.md` for the full version history. Highlights:

- **v2.2.0** — Renamed from `quickclocks` to `clocker`; reverted imputation
  to use `reference_betas.rds`; methylQC sex caller (verbatim).
- **v2.1.x** — kNN imputation, methylQC-style sex inference, embedded
  coefficients, CRAN-compliant cache, missing-probe report, progress bars.
- **v1.0.0** — Initial single-function interface.
