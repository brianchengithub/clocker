# clocker

**Unified Epigenetic Clock Calculator**

A comprehensive R package for calculating 40+ epigenetic clocks from DNA
methylation data with a single function call. Integrates clocks from
multiple research labs into one interface, accepts IDAT files or beta
matrices, and handles imputation, sex inference, cell deconvolution, and
clock projection automatically.

---

## Table of contents

- [Features](#features)
- [Supported platforms](#supported-platforms)
- [First-time setup](#first-time-setup)
- [Quick start](#quick-start)
- [Function arguments](#function-arguments)
- [Output schema](#output-schema)
- [Verbose vs quiet mode](#verbose-vs-quiet-mode)
- [How it works (high level)](#how-it-works-high-level)
- [Imputation](#imputation)
- [Sex inference](#sex-inference)
- [Embedded clock coefficients](#embedded-clock-coefficients)
- [Caching](#caching)
- [Per-sample diagnostics](#per-sample-diagnostics)
- [Performance](#performance)
- [Working on the package itself](#working-on-the-package-itself)
- [Troubleshooting](#troubleshooting)
- [File structure](#file-structure)
- [Citation](#citation)

---

## Features

- **40+ epigenetic clocks** from across the literature in one place
- **Single function call**: `clocker(input)` does everything
- **Many input formats**: IDAT directory, IDAT file paths, beta matrix,
  or `.rds` / `.qs2` / `.csv` / `.tsv` file
- **EPICv2 native**: probe-name normalization with EPIC/450K position
  matching where possible
- **kNN imputation** with reference-mean zero-shot fallback
- **methylQC-style sex inference** (exact algorithm, verified bit-for-bit)
- **EpiDISH cell deconvolution** (RPC + CP, 7 blood cell types)
- **PC-Clocks** with automatic Age/Female fallback when pheno is missing
- **Embedded coefficients** for the simple clocks — no external R-package
  dependency required at runtime once `data-raw/build_coefficients.R`
  has been run
- **CRAN-compliant cache**: `tools::R_user_dir("clocker", "cache")`
- **Per-sample missing-probe report** to spot samples where any clock's
  coverage is degraded
- **Quiet-mode-respecting**: `verbose = FALSE` produces zero stdout,
  with non-default fallbacks surfaced as standard R warnings

## Supported platforms

| Platform        | Probes      | Status            |
| --------------- | ----------- | ----------------- |
| EPICv2 / EPIC+  | ~930,000    | Fully supported   |
| EPIC v1         | ~865,000    | Fully supported   |
| 450K            | ~485,000    | Fully supported   |
| 27K             | ~27,000     | Supported         |
| MSA             | ~285,000    | Supported         |

---

## First-time setup

Setup has three steps. Run them once on each machine where you want to
use clocker.

### Step 1 — install the dependencies

```r
# From the cloned clocker repo, or after downloading install_dependencies.R:
source("install_dependencies.R")
```

This installs the required Bioconductor packages (`sesame`, `EpiDISH`,
`sesameData`) and GitHub-only packages (`DunedinPACE`, `methylCIPHER`,
`EpiMitClocks`), plus optional helpers (`qs2`, `digest`, `progress`,
`matrixStats`, `curl`) that clocker uses when available and falls back
gracefully when they aren't.

### Step 2 — install clocker

```r
pak::pak("brianchengithub/clocker")
```

### Step 3 — build the embedded clock coefficients (recommended)

This step is the one most people miss. clocker ships *without* the
embedded coefficients file because it's generated from
`methylCIPHER` and `EpiMitClocks` at build time, not at install time.
Without this step, clocker still works but pulls coefficients live from
those packages on every call — slower and more brittle.

```r
# After cloning the repo:
setwd("/path/to/clocker/package")
Rscript data-raw/build_coefficients.R
```

What this does:

- Reads coefficient datasets (`Horvath1_CpGs`, `PhenoAge_CpGs`,
  `dataETOC3.l`, etc.) from your installed `methylCIPHER` and
  `EpiMitClocks` packages.
- Saves them as a single xz-compressed file at
  `data/clock_coefficients.rda` (a few hundred KB).
- That file ships with the package the next time you reinstall.

After this step, the simple clocks (Horvath1, Horvath2, Hannum,
PhenoAge, DNAmTL, Lin, Zhang, Bocklandt, Weidner, Vidal-Bralo,
Garagnani, AdaptAge, CausAge, DamAge, SystemsAge, epiTOC2) compute
without needing methylCIPHER or EpiMitClocks installed. PC-Clocks,
DunedinPACE, and the methylCIPHER `calc*` functions still need their
respective backing packages.

If you skip this step, clocker still runs — it just falls back to live
methylCIPHER / EpiMitClocks lookups. You'll see a slightly slower first
run and a few harmless `data set X not found` notes from the upstream
packages.

### Step 4 — verify your install works

A 30-second smoke test that exercises the imputation, sex inference,
and clock projection paths using real CpG IDs:

```r
library(clocker)
packageVersion("clocker")     # confirm a version was loaded

# Use real Horvath1 CpG IDs so the clocks actually compute:
e <- new.env(parent = emptyenv())
utils::data("Horvath1_CpGs", package = "methylCIPHER", envir = e)
horvath_probes <- as.character(get("Horvath1_CpGs", envir = e)$CpG)
horvath_probes <- horvath_probes[!is.na(horvath_probes) &
                                 horvath_probes != "(Intercept)"]

set.seed(42)
all_probes <- unique(c(horvath_probes,
                       paste0("cg", sprintf("%08d", 1:5000))))
betas <- matrix(rbeta(length(all_probes) * 5, 2, 5),
                nrow = length(all_probes), ncol = 5,
                dimnames = list(all_probes, paste0("S", 1:5)))
betas[sample(length(all_probes), 50), 1] <- NA   # tiny missingness

result <- clocker(betas, verbose = TRUE)

# You should see 5 distinct numeric values, not all NA:
result$Clock_Horvath1
```

If `result$Clock_Horvath1` is 5 numbers (not all NA), your install is
healthy.

> **Note**: a smoke test using only fake `cg00000001`-style probe IDs
> won't compute any clocks (no probe overlaps), so the output will look
> empty even when the install is fine. Always seed real CpG IDs into
> the synthetic matrix to actually exercise the projection paths.

---

## Quick start

```r
library(clocker)

# All clocks from an IDAT directory
results <- clocker("/filepath/idats/")

# Or from a pre-processed beta matrix (CpGs x samples, [0, 1] values)
results <- clocker(beta_matrix)

# With phenotype info — pheno$Age and pheno$Sex are auto-coerced
pheno <- data.frame(
  Age = c(45, 67, 22),
  Sex = c("Female", "M", "f"),       # accepts 0/1, M/F, Male/Female (any case)
  row.names = c("S1", "S2", "S3")
)
results <- clocker(beta_matrix, pheno = pheno)

# Quiet mode: no stdout at all
results <- clocker(beta_matrix, verbose = FALSE)
warnings()  # any pheno-fallback alerts go here under verbose=FALSE
```

### Backward-compat aliases (deprecated but still working)

`qclock()`, `calculate_clocks()`, and `run_epigenetic_clocks()` all call
`clocker()` and emit a deprecation warning. Old scripts continue to work.

## Function arguments

```r
results <- clocker(
  input,                                # required
  pheno                   = NULL,       # optional phenotype data.frame
  n_cores                 = NULL,       # NULL = auto-detect
  reference_path          = NULL,       # path to reference_betas.rds
  knn_k                   = 10,         # neighbours for kNN imputation
  knn_zero_shot_threshold = 0.10,       # > this fraction missing -> zero-shot
  missing_report_path     = NULL,       # CSV path for per-sample coverage
  verbose                 = TRUE        # FALSE = silent, alerts -> warnings
)
```

`pheno` recognized columns:

- `Age` (numeric) — used for age-acceleration calculation and as PC-Clocks input
- `Female` (0 / 1) **OR** `Sex` / `sex` — coerced to integer Female with
  `coerce_female_indicator()`. Accepts: `0`, `1`, `"M"`, `"F"`, `"Male"`,
  `"Female"` in any case.

The pheno frame must have rownames matching `colnames(betas)`, OR a
`sample_id` / `Sample_ID` / `id` / `ID` / `Sample` / `sample` column.
clocker will reorder rows to match the betas — it will hard-error if any
sample is missing a pheno row.

## Output schema

A `data.frame` with one row per sample. Column order:

| Group                | Columns                                                                                    |
| -------------------- | ------------------------------------------------------------------------------------------ |
| Sample ID            | `sample_id`                                                                                |
| Sex inference        | `InferredSex`, `InferredSexFlag`, `InferredSexScale`, `chrX_signal`, `chrY_signal`         |
| First-gen clocks     | `Clock_Horvath1`, `Clock_Horvath2`, `Clock_Hannum`, `Clock_PhenoAge`, `Clock_DNAmTL`, ... |
| PC-Clocks            | `PC_Horvath1`, `PC_Horvath2`, `PC_Hannum`, `PC_PhenoAge`, `PC_GrimAge`, `PC_DNAmTL`        |
| Pace / mitotic       | `DunedinPACE`, `epiTOC2_TNSC`, `EpiCMIT`, `MiAge`, `Replitali`                             |
| Specialty            | `Clock_AdaptAge`, `Clock_CausAge`, `Clock_DamAge`, `Clock_SystemsAge`                      |
| Cell composition     | `CellType_RPC_*`, `CellType_CP_*` (B, NK, CD4T, CD8T, Mono, Neutro, Eosino)                |
| Age acceleration     | `Accel_*` (only when `pheno$Age` was supplied)                                             |

Per-sample diagnostics (not in the output frame) are stashed in
`clocker:::.qc_env` after every run — see [Per-sample diagnostics](#per-sample-diagnostics).

## Verbose vs quiet mode

clocker is designed to be quiet by default in scripts and chatty in
interactive use.

**`verbose = TRUE` (default)** — full progress log:

```
=== clocker: Epigenetic Clock Calculator ===
  Platform detected: EPIC
  Using 9 core(s)
  No missing values to impute
  Computing clocks via direct implementations...
    Loaded 29 coefficient datasets
    Direct clocks: Horvath1, Horvath2, Hannum, PhenoAge, DNAmTL, ...
  Estimating cell composition with EpiDISH...
    EpiDISH RPC: 7 cell types estimated
    EpiDISH CP: 7 cell types estimated
  Inferring sex (methylQC algorithm)...
    Sex clustering: threshold=0.707, F=28, M=11, unclear=0
  Computing PC clocks...
    Using cached PC-Clocks data (2.14 GB)
  PC-Clocks: 'Age' not provided; using Horvath2 estimate as Age
  PC-Clocks: 'Female' not provided; using InferredSex (F=1, M/U=0)
    PC clocks added: PCHorvath1, PCHorvath2, ...
=== Done in 12.9 sec ===
```

**`verbose = FALSE`** — zero stdout output; only fallback alerts emitted as warnings:

```r
result <- clocker(beta_matrix, verbose = FALSE)
# (no output)

warnings()
# 1: PC-Clocks: 'Age' not provided; using Horvath2 estimate as Age
# 2: PC-Clocks: 'Female' not provided; using InferredSex (F=1, M/U=0)
```

Use `verbose = FALSE` for batch scripts and pipelines; you'll still hear
about non-default fallbacks via `warnings()`.

## How it works (high level)

1. **Load + validate.** IDAT files go through SeSAMe (`openSesame(prep="QCDPB")`).
   Beta matrices are sanity-checked; out-of-range values are clipped to
   [0, 1]. Robust M-value detection (median + IQR) auto-converts M-values
   when needed without being fooled by a few outlier probes.
2. **EPICv2 normalization.** Suffix-bearing probe names (`cg00000029_TC11`)
   are mapped to base CpG IDs. When multiple replicates exist, the one
   matching the original EPIC v1 / 450K probe by genomic position +
   design + strand is preferred over averaging.
3. **kNN imputation** with reference-mean zero-shot fallback (see below).
4. **methylQC sex inference** — exact port of methylQC's algorithm.
5. **Clock projection.** All available clocks are computed; PC-Clocks
   data (~2 GB) downloads on first use to the cache.
6. **Cell deconvolution.** EpiDISH RPC + CP against `centDHSbloodDMC.m`.
7. **Output assembly.** Stable column ordering with `Clock_*`, `PC_*`,
   `Accel_*`, `CellType_*` prefixes.

## Imputation

Two-pass hybrid:

- **Within-batch kNN.** For each sample with less than
  `knn_zero_shot_threshold` missing-probe fraction (default 10%), the
  `k` most similar OTHER samples in the input batch are found by
  `1 - Pearson` distance over probes both samples have non-missing.
  Each missing probe is imputed as a Gaussian-kernel-weighted mean
  across those neighbours.
- **Zero-shot reference fallback.** Samples above the threshold (or any
  residual probes the kNN neighbours also lack, or single-sample runs)
  fall through to the per-probe means in `reference_betas.rds`.

### `reference_betas.rds`

A **named numeric vector** of mean beta values, keyed by CpG ID. clocker
ships a default in `inst/extdata/reference_betas.rds` (~866K probes
covering the EPIC manifest). For most users, no setup is needed.

To use your own reference (e.g., tissue-specific probe means):

```r
# Build a named numeric vector keyed by CpG ID:
my_ref <- rowMeans(my_reference_beta_matrix, na.rm = TRUE)
saveRDS(my_ref, "/path/to/my_reference_betas.rds")

# Tell clocker to use it (pick one):
options(clocker.reference_betas = "/path/to/my_reference_betas.rds")
Sys.setenv(CLOCKER_REFERENCE_BETAS = "/path/to/my_reference_betas.rds")
# Or pass per-call:
result <- clocker(betas, reference_path = "/path/to/my_reference_betas.rds")
```

Resolution order (first found wins):

1. `reference_path` argument to `clocker()`
2. `options(clocker.reference_betas)`
3. `Sys.getenv("CLOCKER_REFERENCE_BETAS")`
4. `system.file("extdata", "reference_betas.rds", package = "clocker")` (the bundled default)

Data-frame and 1-column-matrix forms are accepted for backward compat,
but the named-numeric-vector form is the canonical one.

## Sex inference

Exact port of methylQC's algorithm:

1. **Per-sample sex signal.** When SigDFs are available (IDAT input),
   total signal intensity (`MG + MR + UG + UR`) is computed over a
   curated panel of 314 chrY probes (PAR-excluded, cross-hyb removed)
   and 3,433 chrX X-inactivation probes (PAR-excluded). For beta-only
   input, median beta over the same probes substitutes.
2. **Threshold optimization.** Sweep candidate chrY thresholds in
   `quantile(chrY, seq(0.15, 0.85, 0.01))` and pick the one minimizing
   total absolute residuals from per-cluster `lm(chrY ~ chrX)` fits.
3. **Confidence bands.** Refit each cluster, compute orthogonal distance
   from each sample to each line, set band thresholds at 5σ.
4. **Classification.** F-band only → F; M-band only → M; both bands →
   tiebreak by chrY threshold; neither band → Unclear.

For batches < 6 samples, embedded reference parameters in
`DEFAULT_SEX_REFERENCE` are used (the data-driven fit needs ≥ 6
samples). Recalibrate for atypical tissues (tumour, placenta, brain) via:

```r
set_sex_reference(
  scale     = "intensity",                # or "beta"
  threshold = 6500,
  male      = list(slope = -0.10, intercept = 7800, sigma = 350),
  female    = list(slope = -0.02, intercept = 1100, sigma = 300)
)
```

## Embedded clock coefficients

clocker bundles clock coefficient tables in `data/clock_coefficients.rda`
to remove runtime dependence on `methylCIPHER` / `EpiMitClocks` for the
simple clocks.

### When to run `data-raw/build_coefficients.R`

- **Once after cloning** the repo (the file isn't committed; it's built
  from your installed methylCIPHER / EpiMitClocks).
- **After updating** methylCIPHER or EpiMitClocks, if you want the
  embedded coefficients to track upstream changes.

### Runtime resolution order

`initialize_clock_coefficients()` tries:

1. The bundled `data/clock_coefficients.rda` (preferred — no external dep)
2. Live `data()` calls into `methylCIPHER` and `EpiMitClocks`
3. Inline R definitions for the smallest clocks (Bocklandt, Garagnani,
   Vidal-Bralo) — always available

After step 1, you can compute Horvath1, Horvath2, Hannum, PhenoAge,
DNAmTL, Lin, Zhang, Zhang2019, Bocklandt, Weidner, VidalBralo,
Garagnani, AdaptAge, CausAge, DamAge, SystemsAge, and epiTOC2 with
**zero external R-package dependencies** beyond base R. PC-Clocks and
DunedinPACE still need their respective packages because they require
training data, not just coefficient tables.

## Caching

clocker caches large files in CRAN-compliant locations. Resolution
order:

1. `getOption("clocker.cache_dir")` if set
2. `Sys.getenv("CLOCKER_CACHE")` if set
3. `tools::R_user_dir("clocker", "cache")` (XDG-compliant on Linux)
4. `tempdir()/clocker` as a final fallback

Cached items:

- `manifests/EPIC.hg38.manifest.qs2` (and platform variants, ~70 MB total)
- `pcclocks/PCClocks_data.qs2` (~2 GB; SHA-256 verified when set)

To find your cache directory:

```r
tools::R_user_dir("clocker", "cache")
```

To redirect the cache (e.g., to a shared scratch volume):

```r
options(clocker.cache_dir = "/scratch/me/clocker_cache")
# or set the CLOCKER_CACHE environment variable in ~/.Renviron
```

## Per-sample diagnostics

After every `clocker()` run, three diagnostic frames live in the package
private environment:

```r
# Imputation: which mode was used per sample, how many neighbours, etc.
clocker:::.qc_env$imputation_info

# Per-clock per-sample probe coverage (also written to disk if you pass
# missing_report_path = "..." to clocker())
clocker:::.qc_env$coverage_log

# Sex caller's per-sample chrX/chrY summary used by the methylQC algorithm
clocker:::.qc_env$sex_signals
```

The imputation frame is the most useful for spotting samples where
imputation may have been unreliable:

```r
diag <- clocker:::.qc_env$imputation_info
diag[diag$mode != "none", ]
# Shows samples where imputation actually fired, with mode, n_neighbors,
# mean_neighbor_dist, n_imputed_knn, n_imputed_zeroshot per sample.
```

## Performance

| Stage                 | Typical cost                           |
| --------------------- | -------------------------------------- |
| First clocker() call  | ~30 sec (downloads PC-Clocks data ~2 GB) |
| Subsequent calls      | ~12-20 sec for 40-60 EPIC samples      |
| First IDAT processing | Slow first time (sesameData cache populates) |
| Memory                | ~5 GB per 100 EPICv2 samples (PC-Clocks projection is the heavy phase) |
| `n_cores`             | Auto-detect uses available RAM and dataset size |

To force a specific core count:

```r
result <- clocker(betas, n_cores = 4)
```

---

## Working on the package itself

If you're modifying clocker (not just using it), this section saves you
some pain.

### The local-vs-GitHub gotcha

`devtools::install()` installs from your local source. `pak::pak()` and
`devtools::install_github()` install from **GitHub**. If you edit
locally and then `pak::pak()` to test, you'll keep getting the
already-pushed version — your local edits won't appear until you commit
and push them.

**Always commit and push before reinstalling from GitHub.**

### Standard edit-test cycle

```r
setwd("/path/to/clocker/")

# 1) Make changes to R/*.R files

# 2) Regenerate documentation from roxygen comments
devtools::document()

# Check for roxygen warnings about missing exports — if you see
# "Objects listed as exports, but not present in namespace", your
# R/ directory is missing source files or has parse errors. Fix that
# before pushing.

# 3) Install locally to test
devtools::install()

# 4) Reload in your session
detach("package:clocker", unload = TRUE)
library(clocker)

# 5) Once working, commit and push (via GitHub Desktop or git CLI)

# 6) Reinstall from GitHub on the test/production machine:
detach("package:clocker", unload = TRUE)
pak::cache_clean()
remove.packages("clocker")
pak::pak("brianchengithub/clocker")
```

### Pre-commit sanity check

This catches the most common "I broke the install" mistakes:

```r
setwd("/path/to/clocker/")

# All R files parse?
for (f in list.files("R", full.names = TRUE)) {
  res <- tryCatch(parse(f), error = function(e) e$message)
  if (!inherits(res, "expression"))
    cat("PARSE ERROR in", f, ":", res, "\n")
}

# Required exports defined?
all_exports <- c("clocker", "qclock", "calculate_clocks",
                 "run_epigenetic_clocks", "set_sex_reference")
sources <- unlist(lapply(list.files("R", full.names = TRUE),
                          readLines, warn = FALSE))
for (fn in all_exports) {
  found <- any(grepl(paste0("^", fn, "\\s*<-\\s*function"), sources))
  cat(sprintf("  %-25s %s\n", fn, if (found) "OK" else "MISSING"))
}

# NAMESPACE has the right exports?
ns <- readLines("NAMESPACE")
n_exports <- length(grep("^export\\(", ns))
cat("NAMESPACE has", n_exports, "exports\n")
```

If any check fails, don't commit — fix it first.

### Bumping the version

`DESCRIPTION` is the single source of truth for the package version.
The startup message (`zzz.R`'s `.onAttach`) reads from it via
`utils::packageVersion()`, so you only need to change one place:

```r
# In R, from the package root:
desc <- readLines("DESCRIPTION")
desc[grepl("^Version:", desc)] <- "Version: 2.2.6"   # whatever
writeLines(desc, "DESCRIPTION")
```

Then commit, push, and reinstall as usual.

### When `pak::pak()` keeps installing the wrong commit

If `pak::pak("brianchengithub/clocker")` reports the same commit hash
even after you've pushed:

```r
pak::cache_clean()                         # nuke pak's local cache
pak::pkg_remove("clocker", lib = .libPaths()[1])
pak::pak("brianchengithub/clocker")
```

If the commit hash still matches the old one, your push didn't reach
GitHub. Check on github.com that the latest commit shows the change you
expect. GitHub Desktop occasionally has staged-but-not-pushed states
that look "done" but aren't.

---

## Troubleshooting

### `EpiDISH not installed` or similar `not installed` messages

```r
BiocManager::install("EpiDISH")
# or
source("install_dependencies.R")  # re-runs the full setup
```

### PC-Clocks data download hangs or times out

Manually download `PCClocks_data.qs2` from
[Zenodo](https://zenodo.org/records/13952402) and place it at:

```r
file.path(tools::R_user_dir("clocker", "cache"), "pcclocks", "PCClocks_data.qs2")
```

clocker will detect the file on next run via its size + SHA-256 check.

### Low probe-overlap warnings

If clocker reports very few probes overlapping with reference matrices
(EpiDISH, clock CpG lists), check that input rownames are standard
Illumina probe IDs (e.g., `cg00000029`). On EPICv2 data, suffixes are
stripped automatically — if you still see the warning, the input may
have been pre-renamed before reaching clocker.

### Clocks producing identical values for every sample

This usually means **the input probe IDs don't match any clock's CpG
list**. The weighted sum reduces to just the intercept, which is the
same for every sample. Check:

```r
# How many of the input probes look like real CpGs?
sum(grepl("^cg|^ch", rownames(betas)))

# How many overlap with a known clock's CpG list?
e <- new.env(parent = emptyenv())
utils::data("Horvath1_CpGs", package = "methylCIPHER", envir = e)
horvath_probes <- get("Horvath1_CpGs", envir = e)$CpG
sum(rownames(betas) %in% horvath_probes)   # should be 300+ for valid input
```

### Verbose mode showing too much chatter from upstream packages

clocker captures `print()`, `cat()`, and `message()` output from
EpiDISH, methylCIPHER, EpiMitClocks, and DunedinPACE. If something
still leaks, please open an issue with `Rscript diagnose_packages.R`
output attached.

### Diagnostic dump

When in doubt, run:

```bash
Rscript diagnose_packages.R > diag.log 2>&1
```

This probes every backend's installation state, runs a synthetic
end-to-end clocker call, and reports cache + reference-resolution
state. Share `diag.log` when reporting issues.

---

## File structure

```
clocker/
├── DESCRIPTION                          # single source of truth for Version
├── NAMESPACE                            # generated by devtools::document()
├── README.md
├── METHODS.md
├── CHANGES.md
├── install_dependencies.R               # one-time dependency setup
├── diagnose_packages.R                  # troubleshooting utility
├── R/
│   ├── clocker.R                        # main function: clocker() + aliases
│   ├── input_processing.R               # IDAT loading, betas validation
│   ├── imputation.R                     # kNN + zero-shot vs. reference_betas.rds
│   ├── sex_inference.R                  # methylQC algorithm (verbatim)
│   ├── sex_probe_lists.R                # 314 chrY + 3,433 chrX probes
│   ├── manifest.R                       # array manifest download/cache
│   ├── clock_engines.R                  # weighted-sum + epiTOC2 + missing-probe report
│   ├── clock_coefficients.R             # embedded coefficients + live fallbacks
│   ├── clock_computation.R              # orchestrator: routes to each clock backend
│   ├── cell_deconvolution.R             # EpiDISH RPC + CP
│   ├── utils.R                          # logging, validation, cache, F/M coercion
│   └── zzz.R                            # package load hooks + private env
├── data-raw/
│   └── build_coefficients.R             # maintainer-only: builds clock_coefficients.rda
├── data/
│   └── clock_coefficients.rda           # generated by data-raw script
└── inst/
    └── extdata/
        └── reference_betas.rds          # default probe-mean reference
```

---

## Citation

If you use clocker, please cite the original publications.

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

See `CHANGES.md` for the full version history.
