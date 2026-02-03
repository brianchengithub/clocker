# 'quickclocks' (Unified Epigenetic Clock Calculator)

A comprehensive R package for calculating 40+ epigenetic clocks from DNA methylation data. Integrates clocks from multiple research labs into a single, easy-to-use interface with a single function call. Accepts many inputs formats -- e.g., directory of IDAT files or beta matrix. 

## Features

- **40+ Epigenetic Clocks** from multiple research groups
- **Single Function Interface** - one call calculates everything
- **Automatic Input Detection** - works with IDAT files or beta matrices
- **EpiDISH Cell Deconvolution** - RPC and CP methods for 7 blood cell types
- **PC-Clocks** - Principal component-adjusted clocks with automatic phenotype inference
- **Optimized Performance** - parallel processing with automatic core detection
- **Smart Caching** - manifests and large data files cached for fast reloading

## Supported Platforms

| Platform | Probes | Status |
|----------|--------|--------|
| EPICv2 / EPIC+ | ~930,000 | ✅ Fully supported |
| EPIC (v1) | ~865,000 | ✅ Fully supported |
| 450K | ~485,000 | ✅ Fully supported |
| 27K | ~27,000 | ✅ Supported |
| MSA | ~285,000 | ✅ Supported |

## Installation

### Option 1: Install as R Package (Recommended)

```r
# Install from GitHub
devtools::install_github("brianchengithub/quickclocks")
```

### Option 2: Source Directly

```r
# Clone the repo, then:
source("install_dependencies.R")  # Install all dependencies
source("R/calculate_clocks.R")    # Load the main function
```

### Dependencies

The install script handles everything automatically:

```r
source("install_dependencies.R")
```

**Key dependencies installed:**
- `sesame` / `sesameData` - IDAT preprocessing
- `EpiDISH` - Cell type deconvolution
- `DunedinPACE` - Pace of aging clock
- `methylCIPHER` - 40+ clocks including PC-Clocks
- `qs2` - Fast data serialization for caching

## Quick Start

```r
# Load the package (or source the file)
library(quickclocks)

# Calculate all clocks from IDAT files
results <- calculate_clocks("/path/to/idat/directory")

# Or from a numeric matrix (of DNA methylation beta values)
results <- calculate_clocks(beta_matrix)

```

### Input Options

**IDAT Files:**
```r
# Point to directory containing IDAT files
results <- calculate_clocks("/data/methylation/idats/")
```
- Supports nested directories
- Automatically finds `*_Grn.idat` and `*_Red.idat` pairs

**Beta Matrix:**
```r
# Provide a pre-processed beta matrix
results <- calculate_clocks(beta_matrix)
```
- Rows: CpG probe IDs (e.g., "cg00000029")
- Columns: Sample IDs
- Values: Beta values [0, 1]

### Options

```r
results <- calculate_clocks(
  input,
  n_cores = NULL,    # NULL = auto-detect optimal cores
  verbose = TRUE     # Print progress messages
)
```

## Output

Returns a `data.frame` with one row per sample and columns for:

| Category | Columns | Description |
|----------|---------|-------------|
| Sample ID | `Sample_ID` | Sample identifiers |
| Cell Types (RPC) | `CellType_RPC_*` | EpiDISH RPC estimates (7 types) |
| Cell Types (CP) | `CellType_CP_*` | EpiDISH CP estimates (7 types) |
| Sex | `InferredSex`, `InferredSex_Numeric` | Chromosome-based inference |
| Age Clocks | `Horvath1`, `Horvath2`, `Hannum`, etc. | Epigenetic age estimates |
| PC Clocks | `PCHorvath1`, `PCHorvath2`, etc. | PC-adjusted clocks |
| Pace | `DunedinPACE` | Pace of biological aging |
| Mitotic | `epiTOC2_TNSC` | Stem cell divisions |

### Example Output

```r
> colnames(results)
 [1] "Sample_ID"           "Horvath1"            "Horvath2"           
 [4] "Hannum"              "PhenoAge"            "DNAmTL"             
 [7] "Lin"                 "Zhang"               "Zhang2019"          
[10] "VidalBralo"          "epiTOC2_TNSC"        "CellType_RPC_B"     
[13] "CellType_RPC_CD4T"   "CellType_RPC_CD8T"   "CellType_RPC_Mono"  
[16] "CellType_RPC_Neu"    "CellType_RPC_NK"     "CellType_RPC_Eos"   
[19] "CellType_CP_B"       "CellType_CP_CD4T"    "CellType_CP_CD8T"   
[22] "CellType_CP_Mono"    "CellType_CP_Neu"     "CellType_CP_NK"     
[25] "CellType_CP_Eos"     "InferredSex_Numeric" "InferredSex"        
[28] "DunedinPACE"         "PCHorvath1"          "PCHorvath2"         
[31] "PCHannum"            "PCPhenoAge"          "PCDNAmTL"           
[34] "PCGrimAge"           "AdaptAge"            "CausAge"            
[37] "DamAge"
```

## Supported Clocks

### Cell Type Deconvolution (EpiDISH)

Uses the `centDHSbloodDMC.m` reference with two methods:

| Method | Description | Output Prefix |
|--------|-------------|---------------|
| **RPC** | Robust Partial Correlations (recommended) | `CellType_RPC_` |
| **CP** | Constrained Projection (Houseman) | `CellType_CP_` |

**Cell Types:** B cells, CD4+ T, CD8+ T, Monocytes, Neutrophils, NK cells, Eosinophils

### First-Generation Age Clocks

| Clock | Year | Description |
|-------|------|-------------|
| `Horvath1` | 2013 | Pan-tissue clock (353 CpGs) |
| `Horvath2` | 2018 | Skin & blood clock |
| `Hannum` | 2013 | Blood-based clock (71 CpGs) |
| `PhenoAge` | 2018 | Mortality-trained clock |
| `GrimAge` | 2019 | Lifespan predictor |
| `DNAmTL` | 2019 | Telomere length estimator |

### PC-Clocks (Principal Component Adjusted)

| Clock | Description |
|-------|-------------|
| `PCHorvath1` | PC-adjusted Horvath multi-tissue |
| `PCHorvath2` | PC-adjusted Horvath skin+blood |
| `PCHannum` | PC-adjusted Hannum |
| `PCPhenoAge` | PC-adjusted PhenoAge |
| `PCGrimAge` | PC-adjusted GrimAge |
| `PCDNAmTL` | PC-adjusted telomere length |

**Note:** PC clocks require a ~2GB data file downloaded automatically from Zenodo on first use.

### Pace of Aging

| Clock | Description |
|-------|-------------|
| `DunedinPACE` | Rate of biological aging (1.0 = average) |

### Mitotic / Cell Division Clocks

| Clock | Description |
|-------|-------------|
| `epiTOC2_TNSC` | Total number of stem cell divisions |

### Specialized Clocks (via methylCIPHER)

| Category | Clocks |
|----------|--------|
| Causality | `AdaptAge`, `CausAge`, `DamAge` |
| Lifestyle | Alcohol, BMI, Smoking predictors |
| Organ-specific | Blood, Brain, Heart, Kidney, Liver, Lung |
| Biological | Inflammation, Immune, Metabolic, Hormone |

## Technical Details

### IDAT Preprocessing Pipeline

When processing IDAT files, uses SeSAMe's `openSesame()`:

1. **NOOB** - Normal-exponential out-of-band background correction
2. **Dye-bias correction** - Non-linear red/green channel normalization
3. **pOOBAH** - P-value with out-of-band array hybridization

### Imputation Strategy

Missing CpG values are imputed using:

1. **Reference-based imputation** - Median values from gold-standard datasets
2. **Zero-shot imputation** - Adds missing clock probes from reference
3. **Fallback** - Row median for probes without reference data

### Sex Inference

Custom chromosome-based method:
- Compares X vs Y chromosome methylation intensity
- Uses platform-specific manifest (auto-downloaded and cached)
- Output: `F` (Female), `M` (Male), or `U` (Unknown)

### PC Clock Phenotype Inference

PC clocks require Age and Sex inputs. When not provided:
- **Age**: Uses Horvath2 (Skin & Blood clock) estimate
- **Sex**: Uses chromosome-based inference

### Caching

The package caches large files for performance:
- **Manifests**: `~/.epigenetic_clock_cache/` (qs2 format)
- **PC Clock data**: Same directory (~2GB, downloaded from Zenodo)

## Performance

### Automatic Thread Optimization

Automatically determines optimal CPU cores based on:
- Dataset size (samples × probes)
- Available CPU cores
- Available RAM

```r
# Auto-detect (recommended)
results <- calculate_clocks(betas)

# Manual override
results <- calculate_clocks(betas, n_cores = 8)
```

### Memory Requirements

| Platform | Memory per 100 samples |
|----------|----------------------|
| 27K | ~0.5 GB |
| 450K | ~2 GB |
| EPIC | ~4 GB |
| EPICv2 | ~5 GB |

### Typical Runtime

- ~20 seconds for 68 EPIC samples (after initial setup)
- First run downloads manifests and PC clock data (~2GB)

## Troubleshooting

### Common Issues

**"EpiDISH not installed"**
```r
BiocManager::install("EpiDISH")
```

**PC clocks timeout during data download**
```r
# Manually download PCClocks_data.qs2 from Zenodo:
# https://zenodo.org/records/13952402
# Place in: ~/.epigenetic_clock_cache/
```

**Low probe overlap warnings**
- Check that your data uses standard Illumina probe IDs (cg########)
- Verify platform detection is correct

### Diagnostic Script

```r
source("diagnose_packages.R")
```

## File Structure

```
EpigeneticClockCalculator/
├── R/
│   ├── calculate_clocks.R      # Main function (single interface)
│   ├── clock_implementations.R # Individual clock calculations
│   ├── clock_availability.R    # Clock definitions and checking
│   ├── imputation.R           # Missing value handling
│   ├── input_handler.R        # IDAT/matrix detection
│   ├── preprocessing.R        # SeSAMe preprocessing
│   └── utils.R                # Helper functions
├── inst/
│   └── extdata/
│       └── reference_betas.rds # Imputation reference data
├── man/                        # Documentation
├── install_dependencies.R      # Dependency installer
├── diagnose_packages.R         # Troubleshooting script
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
└── README.md
```

## Citation

If you use this package, please cite the original publications:

**Core Methods:**
- **EpiDISH**: Teschendorff AE, et al. (2017) *Genome Biology*
- **SeSAMe**: Zhou W, et al. (2018) *Nucleic Acids Research*
- **DunedinPACE**: Belsky DW, et al. (2022) *eLife*
- **PC-Clocks**: Higgins-Chen AT, et al. (2022) *Nature Aging*
- **epiTOC2**: Teschendorff AE (2020) *Genome Biology*
- **methylCIPHER**: Higgins-Chen AT, et al. (2022)

**Original Clocks:**
- **Horvath**: Horvath S (2013) *Genome Biology*
- **Hannum**: Hannum G, et al. (2013) *Molecular Cell*
- **PhenoAge**: Levine ME, et al. (2018) *Aging*
- **GrimAge**: Lu AT, et al. (2019) *Aging*

## License

MIT License - See [LICENSE](LICENSE) file.

## Contributing

Contributions welcome! Please open an issue or submit a pull request.

## Changelog

### v1.0.0
- Single function interface (`calculate_clocks()`)
- EpiDISH cell deconvolution (RPC + CP methods)
- PC clocks with automatic phenotype inference
- Horvath age transformation fix
- qs2 caching for improved performance
- Chromosome-based sex inference
