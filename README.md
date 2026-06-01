# Replication Package: Forestry Leakage in NSW State Forests

**Paper:** Cho, F. H. T., Rhodes, J. R., Maron, M., Lunney, D., Iftekhar, M. S., McIlwee, A., Williams, B. A. (forthcoming). *Leakage effects of establishing protected areas in multi-use forests*

---

## Overview

This package contains the R code and pre-computed data needed to reproduce all figures and summary tables in the paper. The study estimates the causal effect of New South Wales (NSW) state forest revocations on deforestation in adjacent areas using a Before–After Control–Impact (BACI) design with propensity-score matching.

Two leakage channels are estimated:

| Channel | Treatment group | Description |
|---|---|---|
| **buffer** | Grid cells on private land within 10 km of a revoked state forest | Leakage to private native forestry |
| **neighbour** | Grid cells within neighbouring state forests within 10 km | Leakage to state-managed forests |

---

## How to Replicate the Figures

**Requirements:** R ≥ 4.2, with packages: `tidyverse`, `ggplot2`, `patchwork`, `sf`, `glue`, `ggtext`, `igraph`, `scales`, `INLA`, `arrow`

To install INLA:
```r
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))
```

**Steps:**
1. Set R working directory to this folder (`forestry_leakage_rep/`).
2. Run `scripts/05_replicate_figures.R`.
3. All figures are saved to `plots/`.

```r
setwd("path/to/forestry_leakage_rep")
source("scripts/05_replicate_figures.R")
```

The optional `ggmagnify` package adds a zoom inset to the coefficient plot. Install with `remotes::install_github("hughjonesd/ggmagnify")` if desired; it is not required.

---

## Running Scripts 01–03 on HPC

Scripts 01–03 are designed for a PBS/Torque job array. Before submitting jobs, install the required packages to a user-writable library (HPC clusters typically do not allow system-wide installs):

```r
lib_loc <- "~/r_lib"

packages <- c(
  "units", "codetools", "chk", "backports",  # base dependencies
  "dplyr", "tidyr", "readr", "stringr", "purrr", "glue",
  "sf", "terra", "raster", "fasterize",
  "arrow", "geoarrow",
  "exactextractr",
  "MatchIt", "randomForest",
  "INLA"
)

new.packages <- packages[!(packages %in% installed.packages(lib = lib_loc)[, "Package"])]
if (length(new.packages)) {
  install.packages(new.packages, lib = lib_loc, repos = "https://cloud.r-project.org")
}
# INLA must be installed separately:
install.packages("INLA",
  lib = lib_loc,
  repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable")
)
```

Then add `lib.loc = lib_loc` to `library()` calls, or set `.libPaths(lib_loc)` at the top of each script before sourcing.

Each script reads the home directory from the `home` variable near the top — set this to your HPC home path before submitting (e.g. `home <- "/mnt/hpccs01/home/<your-username>"`).

---

## Directory Structure

```
forestry_leakage_rep/
├── scripts/
│   ├── config.R                   Configuration: data and output paths
│   ├── helpers.R                  Shared utility functions
│   ├── ame.R                      Average Marginal Effect (AME) calculation
│   ├── 01_extract_covariates.R    [Reference] Extract pixel-level covariates (HPC)
│   ├── 02_matching_pairing.R      [Reference] Propensity-score matching (HPC)
│   ├── 03_matching_estimation.R   [Reference] BACI estimation via INLA (HPC)
│   ├── 04_process_results.R       [Reference] Aggregate model outputs
│   └��─ 05_replicate_figures.R     ← START HERE: generate all figures
│
├── data/
│   ├── processed_data_for_plotting.RData    Pre-computed model outputs (103 MB)
│   ├── state_forests/
│   │   └── state_forests_revoked_update.gpkg  Revoked SF polygons with revocation year
│   └── focal_state_forest_gpkg/
│       └── focal_sf_year.gpkg              Focal SF polygons with woody veg area estimates
│
├── ext_data/
│   └── *.xlsx                              ABARES log production data (used in Figure 1)
│
└── plots/                                  Output figures (Figures 1–3 pre-populated; run Script 05 for all six)
```

---

## Data Files

| File | Description | Required for |
|---|---|---|
| `data/processed_data_for_plotting.RData` | All model estimates, AME draws, matched data summaries, and spatial objects needed for plotting | Script 05 |
| `data/state_forests/state_forests_revoked_update.gpkg` | NSW state forest polygons with revocation year, area, and woody vegetation attributes | Scripts 04, 05 |
| `data/focal_state_forest_gpkg/focal_sf_year.gpkg` | Focal state forest polygons aggregated by revocation year, with woody vegetation area | Scripts 04, 05 |
| `ext_data/5 Dashboard data tables…xlsx` | ABARES plantation area and log production data (NSW, 1988–2024) | Script 05 (optional) |

### Data not included (large inputs for HPC pipeline)

The following data are required to run Scripts 01–04 from scratch but are **not included** due to size constraints (~100 GB total). They can be obtained from the sources listed below:

| Data | Source | Used in |
|---|---|---|
| SLATS deforestation rasters (1988–2020) | NSW Department of Planning and Environment — [Statewide Landcover and Trees Study](https://datasets.seed.nsw.gov.au/) | Script 01 |
| Woody vegetation annual layers (1988–2020) | Derived from SLATS; contact authors | Script 01 |
| Covariate rasters (25+ layers: land use, ecological condition, climate, etc.) | Various NSW government spatial data portals | Script 01 |
| NSW LGA boundary shapefile | [ABS](https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files) | Script 01 |
| NPWS protected area boundaries | [NSW National Parks and Wildlife Service](https://datasets.seed.nsw.gov.au/) | Script 02 |
| NSW 10 km grid (generated by Script 01) | Generated; stored as `nsw_grid_points.parquet` (1.2 GB) | Scripts 02, 04 |

---

## Script Descriptions

### `scripts/05_replicate_figures.R` — Main replication script

Loads `data/processed_data_for_plotting.RData` and produces all paper figures:

| Figure | File | Description |
|---|---|---|
| Figure 1 | `plots/total_revocations_plt.png` | Bar chart of revoked woody vegetation area by year with NSW log production overlay |
| Figure 2 | `plots/comb_coef_plts.png` | 4-panel: BACI coefficient, total treatment effect (ha), leakage multiplier (%), leakage per ha |
| Figure 3 | `plots/cum_clearing_year_plt.png` | Cumulative deforestation trajectories for treatment vs. control by revocation year |
| Figure 4 | `plots/focal_sf_plt.png` | Spatial map of leakage estimates for each focal state forest |
| Figure 5 | `plots/comb_coef_sf_plt.png` | Combined dashboard: key results and spatial maps |
| Figure 6 | `plots/cumulative_avg_dist_plt.png` | Cumulative weighted-average treatment effect over time |

Also produces:
- `data/leakage_sf_tidy.gpkg` — Tidy spatial export of leakage estimates with credible intervals

### `scripts/04_process_results.R` — Data aggregation [Reference]

Loads raw INLA model RDS files and matching outputs; computes AME, leakage multipliers, and revocation summaries. Saves everything to `data/processed_data_for_plotting.RData`. Requires full HPC output (~50 GB of model files).

### `scripts/03_matching_estimation.R` — BACI estimation [Reference]

Fits negative-binomial INLA models with BACI design for each revocation year and leakage channel. Also runs 1–10 km distance sensitivity analysis. Requires matched pair Parquet files from Script 02.

### `scripts/02_matching_pairing.R` — Propensity-score matching [Reference]

For each revoked state forest, constructs treatment and control groups and runs 19 MatchIt specifications. Selects the best specification using parallel trends, SMD balance, and proportion-retained criteria. Requires the full covariate Parquet files from Script 01.

### `scripts/01_extract_covariates.R` — Covariate extraction [Reference]

Extracts annual deforestation (SLATS), woody vegetation, protected area coverage, and spatial covariates for each 10 km grid cell in NSW. Designed to run as a PBS job array on HPC (1000 jobs). Requires large raster inputs.

### `scripts/helpers.R` — Shared utilities

Contains functions for loading INLA models, computing BACI statistics, pulling and summarising match data, and plotting.

### `scripts/ame.R` — Average Marginal Effect

Contains `generate_ame_irr()`: draws posterior samples from an INLA negative-binomial model and returns the distribution of the total treatment-induced change in cleared area (hectares).

### `scripts/config.R` — Path configuration

Sets `DATA_DIR` and `PLOT_DIR` relative to the package root.

---

## Methodology

The analysis uses a **Difference-in-Differences / BACI** design:

1. Each revoked state forest is a *natural experiment*: clearing is prohibited inside during the revocation event.
2. Grid cells in the leakage zone (buffer: private land; neighbour: adjacent state forests) form the *treatment group*; comparable unaffected cells form the *control group*, selected via propensity-score matching on ecological, climatic, and socioeconomic covariates.
3. The BACI estimand captures the differential change in clearing rate: `(After_Treatment − Before_Treatment) − (After_Control − Before_Control)`.
4. INLA with a negative-binomial family is used to estimate the BACI coefficient, with `log(total_woody_cells)` as an offset.
5. Average Marginal Effects (AME) are computed via posterior sampling from the INLA model and scaled to hectares.

**Matching criteria** (all three must hold for a specification to be selected):
- **Parallel trends**: Bayesian test of pre-treatment trend divergence using INLA (p ≥ 0.05)
- **SMD balance**: All standardised mean differences < 0.25 after matching
- **Proportion retained**: ≥ 80% of treated units matched

---

## Citation

If you use this code or data, please cite:

> Cho, F. (2026). *Leakage from State Forest Revocations in New South Wales: A Quasi-Experimental Analysis.* [Journal TBC].

---

## Contact

For questions about the code or data, please contact the corresponding author or open an issue on the associated repository.
