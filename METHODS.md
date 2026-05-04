# clocker — Methods and Technical Documentation

**Version 2.2.0**

This document describes the methods used in `clocker`, an R package for
calculating 40+ epigenetic clocks from DNA methylation data. It is a full
specification (not a delta to an earlier document) and supersedes any
prior `METHODS.md` from the `quickclocks` package.

---

## 1. Overview

`clocker` consolidates clock implementations from across the literature
into a single function (`clocker()`) that:

1. Accepts IDAT files or pre-processed beta matrices
2. Runs SeSAMe preprocessing for IDAT input
3. Normalizes EPICv2 probe names where needed
4. Imputes missing values via kNN with a zero-shot reference fallback
5. Infers sex via the methylQC algorithm
6. Estimates blood cell composition via EpiDISH
7. Projects all available clocks (first-generation, PC, mitotic, pace)
8. Optionally computes age-acceleration residuals against chronological age
9. Returns a single tidy `data.frame`

All calculations route through a private package environment (`.qc_env`)
so the user's `.GlobalEnv` is never modified.

---

## 2. Input handling

### 2.1 Supported input formats

| Form                                | Treatment                                       |
| ----------------------------------- | ----------------------------------------------- |
| Beta matrix (CpGs × samples)        | Validated, used directly                        |
| `.rds` file (beta matrix)           | `readRDS()`, then validated                     |
| `.qs2` file (beta matrix)           | `qs2::qs_read()` if `qs2` is installed          |
| `.csv` / `.tsv` file (beta matrix)  | `read.csv` / `read.table`, `row.names = 1`      |
| Single IDAT file path               | SeSAMe with basename-derived sample ID          |
| Vector of IDAT file paths           | SeSAMe streaming                                |
| Directory containing IDAT files     | `sesame::openSesame(dir, prep = "QCDPB")`       |

### 2.2 SeSAMe preprocessing pipeline

For IDAT input, `process_idat_files()` calls
`sesame::openSesame(prep = "QCDPB", func = NULL)`:

1. **Q** — Quality masking (multi-mapping, SNP-overlap, cross-hyb)
2. **C** — Channel correction (NOOB)
3. **D** — Dye-bias correction (non-linear)
4. **P** — pOOBAH detection p-values
5. **B** — BMIQ Type-I/Type-II equalization (where applicable)

After SeSAMe, before `getBetas()` is called, `clocker` computes per-sample
chrX and chrY total signal intensities from the SigDFs (the methylQC
canonical input for sex inference) and caches them in `.qc_env`. SigDFs
are then released to keep memory usage bounded.

### 2.3 Validation of beta matrices

`validate_betas()`:

- Requires CpG IDs in rownames and sample IDs in colnames
- Warns if < 50% of rownames look like CpG probes (`^cg|^ch`)
- **Robust M-value detection**: median + IQR rather than min/max — single
  out-of-range probes in 800k don't trigger spurious M→β conversion
- Out-of-range β values within otherwise-valid β-space data are clipped
  to `[0, 1]` with a console message

---

## 3. EPICv2 normalization

EPICv2 uses suffix-bearing probe names (e.g., `cg00000029_TC11`); most
clock coefficient tables, the EpiDISH reference, and DunedinPACE expect
base CpG IDs. `normalize_epicv2_probes()`:

1. Detects suffixed probes via regex `_[TB]C[0-9]{2}$`
2. For CpGs with multiple replicates, calls `resolve_replicates_by_position()`:
   tries (in order) `Probe_beg + DESIGN + strand` → `Probe_beg + DESIGN` →
   `Probe_beg` → `DESIGN + strand`. The replicate matching the original
   EPIC v1 / 450K probe by genomic position + design + strand is preferred.
3. Replicates that cannot be disambiguated by position are averaged.
4. Single-replicate probes are renamed in place.

Manifests are downloaded from the Zhou-lab Infinium V1 annotation
(GitHub) and cached locally. Sex chromosome probes are excluded from
position-resolution attempts to avoid confounding the sex caller
downstream.

---

## 4. Imputation

`clocker` uses a **two-pass hybrid imputation** scheme.

### 4.1 Pass 1: within-batch kNN

For each sample whose missing-probe fraction is at most
`knn_zero_shot_threshold` (default 0.10):

1. **Distance.** For each candidate neighbour (other samples in the input
   batch), compute distance over only those probes that are non-missing
   in both target and candidate. Default metric: `1 - Pearson correlation`,
   robust to additive/multiplicative batch effects. Euclidean is also
   available.

2. **Minimum-coverage check.** If the target has fewer than 50 non-missing
   probes overlapping any candidate, kNN is bypassed and the sample is
   routed to Pass 2.

3. **Neighbour selection.** The `k` reference samples with smallest
   distance are selected (default `k = 10`).

4. **Gaussian-kernel weights.** With bandwidth equal to the median of
   the `k` neighbour distances:

   $$w_i = \frac{\exp(-(d_i / \mathrm{med}(d))^2)}{\sum_j \exp(-(d_j / \mathrm{med}(d))^2)}$$

5. **Per-probe imputation.**

   $$\hat\beta_{p, s} = \sum_{i=1}^{k} w_i\, \beta_{p, r_i}$$

   with the weight renormalized over neighbours non-missing at probe $p$.

### 4.2 Pass 2: zero-shot reference fallback

For samples whose missingness exceeds the threshold, single-sample runs,
or any residual probes that the kNN neighbours also lack:

$$\hat\beta_{p, s} = \bar\beta_p^{\mathrm{ref}}$$

where $\bar\beta_p^{\mathrm{ref}}$ is the per-probe mean from
`reference_betas.rds`. Probes absent from the reference fall back to 0.5.

### 4.3 Reference file: `reference_betas.rds`

A **named numeric vector** of mean β values keyed by CpG ID. For
backward compatibility, `clocker` also accepts a 1-column data frame, a
1-column matrix, or a list with a `$means` element.

Resolution order (first existing wins):

1. `reference_path` argument to `clocker()`
2. `getOption("clocker.reference_betas")`
3. `Sys.getenv("CLOCKER_REFERENCE_BETAS")`
4. `system.file("extdata", "reference_betas.rds", package = "clocker")`

### 4.4 Diagnostics

After every run, `.qc_env$imputation_info` contains a per-sample
`data.frame` with: `sample_id`, `n_probes_total`, `n_missing_input`,
`pct_missing_input`, `mode` ∈ {`knn`, `zero_shot`, `knn_failed_zeroshot`,
`fallback_constant`, `none`}, `n_neighbors_used`, `mean_neighbor_dist`,
`n_imputed_knn`, `n_imputed_zeroshot`. This is the primary tool for
spotting samples where imputation may be unreliable.

### 4.5 References

- Troyanskaya O et al. (2001) *Bioinformatics* 17:520 (kNN imputation)
- Hechter E et al. (2018) *BMC Bioinformatics* 19:354
- Lent S et al. (2019) *Epigenetics* 14:1

---

## 5. Cell-type deconvolution

`estimate_cell_composition()` uses `EpiDISH::epidish()` against the
`centDHSbloodDMC.m` reference (7 blood cell types: B, NK, CD4T, CD8T,
Mono, Neutro, Eosino). Both methods are run:

- **RPC** (robust partial correlations) — default, recommended
- **CP** (constrained projection) — Houseman-style, also returned

Output columns: `CellType_RPC_*` and `CellType_CP_*`. Probe overlap
with the reference is reported in the log; < 100 overlapping probes
triggers a low-reliability warning.

For non-blood tissues, run EpiDISH manually with an appropriate reference
matrix; cell-type estimation can be skipped without affecting the rest
of `clocker`.

### References

- Teschendorff AE et al. (2017) *BMC Bioinformatics* 18:105 (EpiDISH)
- Houseman EA et al. (2012) *BMC Bioinformatics* 13:86 (CP)
- Reinius LE et al. (2012) *PLOS One* 7:e41361 (blood reference)

---

## 6. Sex inference (methylQC algorithm — exact)

`clocker`'s sex caller is a verbatim port of methylQC's
`compute_sex_intensities_single()` and `plot_sex_check_optimal()` (MIT
license, by Brian Cheng). It has been verified bit-for-bit equivalent to
methylQC on synthetic data.

### 6.1 Stage 1 — per-sample sex signal

Two curated probe sets, embedded in `R/sex_probe_lists.R` (extracted
from `sesameDataGet("EPIC.probeInfo")`):

- `chrY_clean` — 314 non-PAR chrY probes with cross-hybridizing probes removed
- `chrX_xlinked` — 3,433 non-PAR X-inactivation chrX probes

Per sample:

- **IDAT input** (preferred): total signal intensity over each probe set,
  $\mathrm{TI} = \mathrm{MG} + \mathrm{MR} + \mathrm{UG} + \mathrm{UR}$
  with `NA → 0`. Per sample,
  $\mathrm{chrX} = \mathrm{med}_{p \in \mathrm{chrX\_xlinked}}(\mathrm{TI}_p)$
  and analogously for chrY.
- **Beta matrix input**: median β over the same probes, used as a proxy.

Both require ≥ 10 matching probes; otherwise the value is `NA`.

The `InferredSexScale` output column records which mode was used.

### 6.2 Stage 2 — data-driven threshold optimization

Candidate thresholds: $\tau \in \mathrm{quantile}(\mathrm{chrY},\,
\mathrm{seq}(0.15, 0.85, 0.01))$.

For each $\tau$ requiring at least 3 samples in each subgroup:

- Females = samples with chrY ≤ τ; fit $\mathrm{lm}(\mathrm{chrY} \sim \mathrm{chrX})$
- Males = samples with chrY > τ; fit analogously
- $\mathrm{cost}(\tau) = \sum |\mathrm{residuals}_F| + \sum |\mathrm{residuals}_M|$

Best $\tau$ = $\arg\min \mathrm{cost}$. If no valid τ is found,
$\tau$ defaults to `median(chrY)`.

### 6.3 Stage 3 — confidence bands

Refit per-cluster regressions on the final F/M assignment. Orthogonal
distance from a point $(x, y)$ to a fitted line with slope $a$ and
intercept $c$:

$$d(x, y) = \frac{|a x - y + c|}{\sqrt{a^2 + 1}}$$

Cluster-specific $\sigma_F, \sigma_M$ = standard deviations of the
in-cluster orthogonal distances. Band thresholds:
$\mathrm{thr}_F = 5\,\sigma_F$, $\mathrm{thr}_M = 5\,\sigma_M$.

### 6.4 Stage 4 — classification

| $d_F \le \mathrm{thr}_F$ | $d_M \le \mathrm{thr}_M$ | Call                                  | Flag                          |
| ------------------------ | ------------------------ | ------------------------------------- | ----------------------------- |
| Yes                      | No                       | F                                     | `F_band_only`                 |
| No                       | Yes                      | M                                     | `M_band_only`                 |
| Yes                      | Yes                      | tiebreak: `chrY ≤ τ → F` else `M`     | `both_bands_chrY_tiebreak`    |
| No                       | No                       | Unclear                               | `outside_both_bands`          |

### 6.5 Small-batch fallback

methylQC's algorithm requires ≥ 6 input samples for the data-driven fit.
For batches smaller than that (or pathological fits), `clocker` falls
back to embedded reference-cluster regressions in `DEFAULT_SEX_REFERENCE`.
The same Stage 3/4 logic runs with reference parameters. Default
parameters are calibrated for whole-blood EPIC; tissues with atypical
X-inactivation should be recalibrated:

```r
set_sex_reference(
  scale     = "intensity",  # or "beta"
  threshold = 6500,
  male      = list(slope = -0.10, intercept = 7800, sigma = 350),
  female    = list(slope = -0.02, intercept = 1100, sigma = 300)
)
```

### 6.6 Output columns

- `chrX_signal`, `chrY_signal` — raw inputs to the classifier
- `InferredSex` — `F`, `M`, or `Unclear`
- `InferredSexFlag` — diagnostic flag from the table above
- `InferredSexScale` — `intensity` (IDAT) or `beta` (matrix input)

### 6.7 References

- methylQC: https://github.com/brianchengithub/methylQC (MIT, by Brian Cheng)
- McCartney DL et al. (2016) *Genom Data* 9:22 (cross-hyb probe lists)
- Zhou W et al. (2017) *NAR* 45:e22 (Infinium probe annotation)

---

## 7. Phenotype handling

### 7.1 Pheno frame validation and alignment

`validate_and_align_pheno()`:

- Accepts `sample_id` / `Sample_ID` / `id` / `ID` / `Sample` / `sample`
  as ID columns; otherwise uses `rownames(pheno)`.
- Hard-errors if any input sample has no pheno row (preventing the
  silent-correctness bug from earlier versions).
- **Reorders pheno rows to match `colnames(betas)`.**
- Coerces `Female` / `Sex` / `sex` columns via `coerce_female_indicator()`.

### 7.2 Female / Sex coercion

`coerce_female_indicator()` accepts:

- `0` / `1` numeric
- `"M"` / `"F"` / `"Male"` / `"Female"` (any case, with whitespace)

Conversion: `F` / `Female` / `1` → `1L`; `M` / `Male` / `0` → `0L`.
Unknown / `0.5` / `"U"` / empty / `NA` → `NA_integer_` with a warning.

### 7.3 Age fallback for PC-Clocks

PC-Clocks require Age and Female. When either is missing:

- **Age missing** → use `Horvath2` estimate; if that's also unavailable,
  fall back to `Horvath1`. Either fallback emits an explicit
  `message()` so the user knows which proxy is in use.
- **Female missing** → use `InferredSex` (`F → 1`, `M / Unclear → 0`);
  if sex inference also failed, default to `0` (Male) with a console
  message.

The resulting Age and Female values are also reported in the log so the
user can audit them.

---

## 8. Clock projection

### 8.1 Generic weighted-sum engine

`calc_weighted_sum_clock()` handles the common shape of clock
coefficient tables: `intercept + sum(beta * weight)` with optional
output transformation. It auto-detects:

- **CpG column**: `CpG`, `CpGmarker`, `probe`, `Probe`, `cpg`, `ID`,
  `ProbeID`, `probe_id`, `CpG_ID`, `Marker`, or the first character/
  factor column.
- **Weight column**: `Coefficient`, `CoefficientTraining`, `weight`,
  `Weight`, `coef`, `beta`, `Beta`, `Effect`, `effect`, or the first
  numeric non-CpG column.

Coverage tracking: per clock per sample, `n_cpgs_total` and
`n_cpgs_present` are recorded so the missing-probe report
(`write_missing_probe_report()`) can later quantify how complete each
clock's input was for each sample.

If multiple intercept rows exist, a warning is emitted and the first is
used. Probe-coverage attributes (`n_cpgs_total`, `n_cpgs_used_max`) are
attached to the return value.

### 8.2 Horvath transformation

For Horvath1 and Horvath2, the raw weighted sum is mapped to age in
years via the piecewise transform from Horvath (2013):

$$\mathrm{age}(x) = \begin{cases} (a + 1) \exp(x) - 1 & x < 0 \\ (a + 1) x + a & x \ge 0 \end{cases}$$

with adult age $a = 20$.

### 8.3 First-generation clocks

| Clock         | Source                 | Output         |
| ------------- | ---------------------- | -------------- |
| Horvath1      | methylCIPHER / embedded | `Clock_Horvath1` |
| Horvath2      | methylCIPHER / embedded | `Clock_Horvath2` |
| Hannum        | methylCIPHER / embedded | `Clock_Hannum` |
| PhenoAge      | methylCIPHER / embedded | `Clock_PhenoAge` |
| DNAmTL        | methylCIPHER / embedded | `Clock_DNAmTL` |
| Lin           | methylCIPHER / embedded | `Clock_Lin` |
| Zhang         | methylCIPHER / embedded | `Clock_Zhang` |
| Zhang2019     | methylCIPHER / embedded | `Clock_Zhang2019` |
| Bocklandt     | methylCIPHER / embedded / inline fallback | `Clock_Bocklandt` |
| Weidner       | methylCIPHER / embedded | `Clock_Weidner` |
| VidalBralo    | methylCIPHER / embedded / inline fallback | `Clock_VidalBralo` |
| Garagnani     | methylCIPHER / embedded / inline fallback | `Clock_Garagnani` |
| AdaptAge      | methylCIPHER            | `Clock_AdaptAge` |
| CausAge       | methylCIPHER            | `Clock_CausAge` |
| DamAge        | methylCIPHER            | `Clock_DamAge` |
| SystemsAge    | methylCIPHER            | `Clock_SystemsAge` |

The `inline fallback` clocks (Bocklandt, Garagnani, VidalBralo) have
hard-coded coefficient tables in `clock_coefficients.R` so they work
even when both methylCIPHER and the embedded `.rda` are unavailable.

### 8.4 PC-Clocks

`compute_pc_clocks()` calls `methylCIPHER::calcPCClocks()` after:

1. Loading `PCClocks_CpGs` into `.qc_env`
2. Resolving Age and Female (with Horvath fallback alerts as in §7.3)
3. Verifying / downloading PC-Clocks training data (~2 GB) from Zenodo
   to `tools::R_user_dir("clocker", "cache")/pcclocks/PCClocks_data.qs2`
4. SHA-256 verification (when `PCCLOCKS_SHA256` is set) plus a minimum
   file-size check (≥ 1 GB)

PC-Clocks output: `PC_Horvath1`, `PC_Horvath2`, `PC_Hannum`,
`PC_PhenoAge`, `PC_GrimAge`, `PC_DNAmTL`.

When methylCIPHER expects coefficient datasets in `.GlobalEnv`, `clocker`
mirrors them temporarily and removes them via `on.exit()` so the user's
workspace stays clean.

### 8.5 DunedinPACE

`compute_dunedin_pace()` calls `DunedinPACE::PACEProjector(betas)`.
Single-sample inputs are temporarily duplicated (`PACEProjector` requires
≥ 2 samples) and the duplicate is dropped from the output. Output:
`DunedinPACE`.

### 8.6 Mitotic clocks

`compute_mitotic_clocks()` invokes `EpiMitClocks::EpiMitClocks(data.m =
betas)` (or `epiTOC2`, etc.) after loading the required datasets
(`dataETOC3`, `cugpmitclockCpG`, `epiTOCcpgs3`, `estETOC3`,
`EpiCMITcpgs`, `Replitali`) into `.qc_env` and mirroring to `.GlobalEnv`
under `on.exit()`.

`epiTOC2` is also computed directly from coefficients in
`calc_epitoc2_direct()` if `dataETOC3.l` is available. Output:
`epiTOC2_TNSC`, `EpiCMIT`, `MiAge`, `Replitali`.

---

## 9. Embedded clock coefficients

To minimize external dependencies, coefficient tables for the simple
clocks are bundled with the package.

### 9.1 Build step (maintainer-only)

```bash
Rscript data-raw/build_coefficients.R
```

This script extracts datasets from installed `methylCIPHER` and
`EpiMitClocks` packages and saves them to `data/clock_coefficients.rda`
(xz-compressed, ~ a few hundred KB).

### 9.2 Runtime resolution order

`initialize_clock_coefficients()` tries:

1. The bundled `data/clock_coefficients.rda`
2. Live `data()` calls into `methylCIPHER` and `EpiMitClocks` (if installed)
3. Inline R definitions for Bocklandt, Garagnani, VidalBralo

After step 1 has been performed once by the maintainer, end users can
compute Horvath1, Horvath2, Hannum, PhenoAge, DNAmTL, Lin, Zhang,
Zhang2019, Bocklandt, Weidner, VidalBralo, Garagnani, AdaptAge, CausAge,
DamAge, SystemsAge, and epiTOC2 with **zero external R-package
dependencies** (beyond base R and EpiDISH for cell deconvolution).

PC-Clocks training data (~2 GB) cannot be embedded due to package size
limits and remains an on-demand Zenodo download (see §10).

---

## 10. Caching

All cached files live in CRAN-compliant locations resolved by:

1. `getOption("clocker.cache_dir")` if set
2. `Sys.getenv("CLOCKER_CACHE")` if set
3. `tools::R_user_dir("clocker", "cache")` (XDG-compliant on Linux)
4. `tempdir()/clocker` as a final fallback

Cached items:

- `manifests/<PLATFORM>.hg38.manifest.qs2` — Zhou-lab Infinium V1
  annotation, downloaded on demand
- `pcclocks/PCClocks_data.qs2` — PC-Clocks training data, ~2 GB, SHA-256
  verified when `PCCLOCKS_SHA256` is set in `R/clock_computation.R`

Corrupt or undersized files are deleted and re-downloaded with
exponential-backoff retry.

---

## 11. Output schema

The `data.frame` returned by `clocker()` has stable column ordering:

1. `sample_id`
2. Sex inference: `InferredSex`, `InferredSexFlag`, `InferredSexScale`,
   `chrX_signal`, `chrY_signal`
3. Cell composition: `CellType_RPC_*`, `CellType_CP_*`
4. First-generation clocks: `Clock_Horvath1`, `Clock_Horvath2`,
   `Clock_Hannum`, `Clock_PhenoAge`, ...
5. PC-Clocks: `PC_Horvath1`, `PC_Horvath2`, `PC_Hannum`, `PC_PhenoAge`,
   `PC_GrimAge`, `PC_DNAmTL`
6. Pace and mitotic: `DunedinPACE`, `epiTOC2_TNSC`, `EpiCMIT`, `MiAge`,
   `Replitali`
7. Specialty: `Clock_AdaptAge`, `Clock_CausAge`, `Clock_DamAge`,
   `Clock_SystemsAge`
8. Age acceleration: `Accel_*` (only when `pheno$Age` is available)

Per-sample diagnostics (not in the output frame) are stashed in
`.qc_env`:

- `imputation_info` — see §4.4
- `coverage_log` — per-clock per-sample probe coverage (used by the
  optional `missing_report_path` CSV writer)
- `sex_signals` — per-sample chrX/chrY summary used by the sex caller

---

## 12. Resource estimation

`determine_optimal_cores()` picks an effective `n_cores` based on:

- Available physical cores from `parallel::detectCores(logical = FALSE)`
- Available RAM from `/proc/meminfo` (Linux), `sysctl hw.memsize` (macOS),
  or `wmic` (Windows), with a fallback of 8 GB if probing fails
- Per-sample memory estimate of $n_{\mathrm{probes}} \times 8 \times 3$
  bytes (input + intermediate + scratch). The `× 3` factor accounts for
  the input matrix, an intermediate copy used during clock projection,
  and scratch space for matrix-multiply temporaries.

The user can override with `n_cores = N`; the function emits a warning
if the user request exceeds available resources.

---

## 13. Logging and diagnostics

`log_msg()` is the single console-output helper. All progress messages
go through it and respect the `verbose` argument.

Per-sample loops (kNN imputation, clock projection) use a unified
progress bar (`make_progress_bar()`) that prefers `progress::progress_bar`
when installed and falls back to `utils::txtProgressBar` otherwise.

Per-sample, per-clock missing-probe coverage can be written to CSV via
`missing_report_path = "..."`. Two files are produced: a long-form
(`sample_id`, `clock`, `n_cpgs_total`, `n_cpgs_present`, `pct_missing`)
and a wide-form companion (samples × clocks of `pct_missing`).

---

## 14. Limitations and caveats

- **PC-Clocks training data is large** (~2 GB) and requires an internet
  connection on first use.
- **Cell deconvolution defaults to a blood reference.** Non-blood tissues
  produce numerically-valid but biologically-meaningless `CellType_*`
  fractions; users should disable or rerun with a tissue-appropriate
  reference.
- **methylQC sex caller defaults are blood-tuned** for the small-batch
  fallback. Recalibrate via `set_sex_reference()` for tumor, placenta,
  or brain.
- **Imputation accuracy depends on `reference_betas.rds`.** A reference
  built from heavily blood-skewed data will impute non-blood inputs less
  accurately. The diagnostic frame in `.qc_env$imputation_info` is the
  canary for this.
- **Cross-platform comparability** of clocks is best when input is EPIC
  or 450K (training-set platforms for nearly all clocks). EPICv2 is
  supported via probe-name normalization; comparability is maintained
  for the resolved (position-matched) probes but not for those that had
  to be averaged across replicates.

---

## 15. References (consolidated)

**Preprocessing**

- Zhou W et al. (2018) *NAR* 46:e123 (SeSAMe)
- Triche T et al. (2013) *NAR* 41:e90 (NOOB)
- Teschendorff AE et al. (2013) *Bioinformatics* 29:189 (BMIQ)

**Cell deconvolution**

- Houseman EA et al. (2012) *BMC Bioinformatics* 13:86
- Reinius LE et al. (2012) *PLOS One* 7:e41361
- Teschendorff AE et al. (2017) *BMC Bioinformatics* 18:105

**Sex inference**

- methylQC: https://github.com/brianchengithub/methylQC (Brian Cheng, MIT)
- McCartney DL et al. (2016) *Genom Data* 9:22

**Clocks**

- Horvath S (2013) *Genome Biology* 14:R115
- Horvath S et al. (2018) *Aging* 10:1758 (Skin & Blood)
- Hannum G et al. (2013) *Molecular Cell* 49:359
- Levine ME et al. (2018) *Aging* 10:573 (PhenoAge)
- Lu AT et al. (2019) *Aging* 11:303 (GrimAge)
- Lu AT et al. (2019) *Aging* 11:5895 (DNAmTL)
- Belsky DW et al. (2022) *eLife* 11:e73420 (DunedinPACE)
- Higgins-Chen AT et al. (2022) *Nature Aging* 2:644 (PC-Clocks)
- Teschendorff AE (2020) *Genome Medicine* 12:56 (epiTOC2)
- Lin Q et al. (2016) *Aging* 8:394
- Zhang Y et al. (2017) *Genome Medicine* 9:21
- Zhang Q et al. (2019) *Genome Medicine* 11:54
- Bocklandt S et al. (2011) *PLOS One* 6:e14821
- Weidner CI et al. (2014) *Genome Biology* 15:R24
- Vidal-Bralo L et al. (2016) *Frontiers in Genetics* 7:126
- Garagnani P et al. (2012) *Aging Cell* 11:1132
- Ying K et al. (2024) *Nature Aging* 4:231 (AdaptAge / CausAge / DamAge)

**Imputation**

- Troyanskaya O et al. (2001) *Bioinformatics* 17:520
- Hechter E et al. (2018) *BMC Bioinformatics* 19:354
- Lent S et al. (2019) *Epigenetics* 14:1

**EPICv2**

- Kaur D et al. (2023) *Genome Biology* 24:101
