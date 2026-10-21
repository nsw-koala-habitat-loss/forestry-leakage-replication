# Replication Package: Forestry Leakage in NSW State Forests

**Paper:** Cho, F. H. T., Rhodes, J. R., Maron, M., Lunney, D., Iftekhar, M. S., McIlwee, A., Williams, B. A. (forthcoming). *Leakage effects of establishing protected areas in multi-use forests*

---

## Overview

This package contains the R code and pre-computed data needed to reproduce all figures and summary tables in the paper. The study estimates the causal effect of New South Wales (NSW) state forest revocations on deforestation in adjacent areas using a Before–After Control–Impact (BACI) design with propensity-score matching.

Two leakage channels are estimated:

| Channel | Treatment group | Description |
|---|---|---|
| **buffer** | Grid cells on private land within 10 km of a revoked state forest | Leakage to private land |
| **neighbour** | Grid cells within neighbouring state forests within 10 km | Leakage to other state forests |

---

## How to Replicate the Figures

**Requirements:** R ≥ 4.2, with packages: `tidyverse`, `ggplot2`, `patchwork`, `sf`, `glue`, `ggtext`, `igraph`, `scales`, `INLA`, `arrow`. Script 06b additionally needs `gt` and `openxlsx` (for the formatted `.xlsx`; `writexl` gives an unformatted fallback); `ggmagnify` is optional (magnifier insets in Figures 2/5).

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
│   ├── 02b_woody_vegetation_areas.R [Reference] Woody-veg areas + on-site avoided clearing (HPC)
│   ├── 03_matching_estimation.R   [Reference] BACI estimation via INLA (HPC)
│   ├── 04_process_results.R       [Reference] Aggregate model outputs
│   ├── 05_replicate_figures.R     ← START HERE: generate all figures
│   ├── 06_figure5_redesign.R      Alternative 2-panel Figure 5
│   ├── 06b_figure5_derivation_table.R  Leakage-multiplier derivation table + heatmap
│   ├── 06c_figure5_prototype_absolute.R Prototype absolute-ha panels for Figure 5
│   ├── 07_figure4_parallel_trends.R    Re-based Figure 4 parallel-trends check (R&R)
│   └── 08a_export_match_map_inputs.R   [Reference] Export inputs for the matching-design map (authors only)
│
├── data/
│   ├── processed_data_for_plotting.RData    Pre-computed model outputs (103 MB)
│   ├── nsw_grid_points.parquet              NSW 10 km grid cell centroids (1.2 GB)
│   ├── nsw_grid_points_README.md            Column reference for nsw_grid_points.parquet
│   ├── state_forests/
│   │   └── state_forests_revoked_update.gpkg  Revoked SF polygons with revocation year
│   ├── focal_state_forest_gpkg/
│   │   └── focal_sf_year.gpkg              Focal SF polygons: woody veg area + on-site avoided clearing (from Script 02b)
│   └── matching_data/
│       ├── README.md                        Column reference for matching data files
│       └── parquet/
│           └── matching_data_pixel_{1..1000}.parquet  Pixel-level covariate panel (1,000 files)
│
├── ext_data/
│   └── *.xlsx                              ABARES log production data (used in Figure 1)
│
└── plots/                                  Output figures (Figures 1–3 pre-populated; run Script 05 for all)
```

---

## Data Files

| File | Description | Required for |
|---|---|---|
| `data/processed_data_for_plotting.RData` | All model estimates, AME draws, matched data summaries, and spatial objects needed for plotting. See [Contents of processed_data_for_plotting.RData](#contents-of-processed_data_for_plottingrdata) below for a full inventory. | Script 05 |
| `data/state_forests/state_forests_revoked_update.gpkg` | NSW state forest polygons with revocation year, area, and woody vegetation attributes. Source: [NSW Forestry Corporation — State Forests of NSW](https://data-fcnsw.opendata.arcgis.com/datasets/bbf614f12b58432a85ab680162a0c9e0_0/explore?location=-32.844559%2C148.895110%2C7) (publicly available; [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/)) | Scripts 04, 05 |
| `data/focal_state_forest_gpkg/focal_sf_year.gpkg` | Focal state forest polygons aggregated by revocation year, with woody vegetation area | Scripts 04, 05 |
| `data/nsw_grid_points.parquet` | NSW 10 km regular grid cell centroids (712,714 cells; 1.2 GB) generated by Script 01. Each row is a 150 m buffered polygon around a grid centroid. See `data/nsw_grid_points_README.md` for details. | Scripts 02, 04 |
| `data/matching_data/parquet/matching_data_pixel_{1..1000}.parquet` | Pixel-level covariate and deforestation panel for all NSW 10 km grid cells, split into 1,000 files (one per HPC job). See `data/matching_data/README.md` for full column reference. | Script 02 |
| `ext_data/5 Dashboard data tables…xlsx` | ABARES plantation area and log production data (NSW, 1988–2024) | Script 05 (optional) |

### Data not included (large inputs for HPC pipeline)

The following raw input data are required to run Scripts 01–04 from scratch but are **not included** due to size constraints (~100 GB total). They can be obtained from the sources listed below:

| Data | Source | Used in |
|---|---|---|
| SLATS deforestation rasters (1988–2014) | NSW Department of Climate Change, Energy, the Environment and Water — [Statewide Landcover and Trees Study](https://datasets.seed.nsw.gov.au/) ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); © State Government of NSW and NSW Department of Climate Change, Energy, the Environment and Water) | Script 01 |
| Woody vegetation annual layers (1988–2014) | Derived from SLATS; contact authors | Script 01 |
| Covariate rasters (25+ layers: land use, ecological condition, climate, etc.) | Various NSW government spatial data portals | Script 01 |
| NSW LGA boundary shapefile | [ABS](https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files) | Script 01 |
| NPWS protected area boundaries | [NSW National Parks and Wildlife Service](https://datasets.seed.nsw.gov.au/) | Script 02 |

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
| Figure 6 | `plots/cumulative_avg_dist_plt.png` | Cumulative and distribution of the weighted-average treatment effect per ha |
| Figure 6 (companion) | `plots/cumulative_avg_leakage_plt.png` | Distribution of the average leakage multiplier (leakage relative to on-site avoided clearing) |

Also produces:
- `data/leakage_sf_tidy.gpkg` — Tidy spatial export of leakage estimates with credible intervals

### `scripts/06_figure5_redesign.R` — Alternative Figure 5

Runs from the same `data/processed_data_for_plotting.RData` as Script 05. Produces a 2-panel alternative to the Figure 5 dashboard — total leakage effect (ha) and leakage multiplier (%), with 95%-CI significance shown by point fill (`plots/figure5_redesign.png`). Reuses Script 04's `leakage_multiplier_df`, so the numbers match the main Figure 5.

### `scripts/06b_figure5_derivation_table.R` — Leakage-multiplier derivation table

Companion to Script 06. Writes a step-by-step derivation table (`plots/figure5_derivation_table.html` / `.png`, plus a formatted `.xlsx` via `openxlsx`, or a plain sheet via `writexl`) and a year × channel heatmap (`plots/figure5_derivation_heatmap.png`). Runs from `data/processed_data_for_plotting.RData`; needs `gt`.

### `scripts/07_figure4_parallel_trends.R` — Re-based Figure 4 (R&R)

Draws a re-based version of Figure 4: both arms set to 0 the year before revocation, y-axis cumulative clearing as a % of plot woody vegetation — the quantity the parallel-trends test evaluates. Panels are flagged where the selected matched specification fails the pre-trend criterion. Runs from `data/processed_data_for_plotting.RData` (`plots/fig4_parallel_trends.png`).

### `scripts/06c_figure5_prototype_absolute.R` — Absolute-ha Figure 5 prototypes

Candidate replacements for panel A of Script 06 that keep the leakage effect in hectares (cumulative and per-year variants). Runs from the shipped `.RData`; writes `plots/figure5_absolute_prototypes.png`.

### `scripts/08a_export_match_map_inputs.R` — Matching-design map inputs [Reference]

Run once by the authors on the full analysis outputs; it uses hard-coded HPC paths and cannot be run from this package. It exports the cropped layers for the matching-design map to `data/match_map/`. The companion plotting script (Script 08) and `data/match_map/` are not yet included in this package.

### `scripts/04_process_results.R` — Data aggregation [Reference]

Loads raw INLA model RDS files and matching outputs; computes AME, leakage multipliers, and revocation summaries. Saves everything to `data/processed_data_for_plotting.RData`. Requires full HPC output (~50 GB of model files).

### `scripts/03_matching_estimation.R` — BACI estimation [Reference]

Fits negative-binomial INLA models with BACI design for each revocation year and leakage channel. Also runs 1–10 km distance sensitivity analysis and writes cumulative-clearing trajectories for all 19 matching specs (Figure 3 input). Requires matched pair Parquet files from Script 02 and `focal_sf_year.gpkg` from Script 02b.

### `scripts/02b_woody_vegetation_areas.R` — Woody-vegetation areas [Reference]

Ported from the human-written `woody_veg_calc.R`. From the annual woody-vegetation and SLATS rasters, computes for each focal state forest: woody vegetation inside the forest at the revocation year (`woody_area_start_m2`), on-site avoided clearing (`five_year_clearing_ha` = class-4 clearing in the 5 years before revocation minus the 5 years after), and woody vegetation inside the 10 km leakage belt (`year_run_df.csv`). Writes `focal_sf_year.gpkg` (shipped pre-computed). Requires the full upstream rasters.

### `scripts/02_matching_pairing.R` — Propensity-score matching [Reference]

For each revoked state forest, constructs treatment and control groups and runs 19 MatchIt specifications. Selects the best specification using parallel trends, SMD balance, and proportion-retained criteria. Requires the full covariate Parquet files from Script 01.

> **Protected-area screen (updated 2026-08).** The PA screen now uses the layer covering the end of the BACI panel — `PA_{revocation_year + 6}` (the annual layers are cumulative: `PA_{y}` = gazetted on/before `y − 1`), was `+ 4` — and is applied to the **neighbour** channel as well as the buffer channel; the neighbour donor-pool "intact state forest" cutoff moved from `revocation_year + 4` to `+ 5` to match. This changes no leakage-multiplier median by more than ~1.1 pp and flips no significance call (buffer 0.6 % / neighbour 1.5 % of matched plots dropped). **`data/processed_data_for_plotting.RData` and `plots/` reflect this correction**, applied post-hoc — protected plots were dropped from the shipped matched pairs and the BACI + AME models re-fit — since re-running the matching itself needs the HPC raster/shapefile inputs. `data/pa_correction_*.csv` records the dropped-plot counts and the shipped-vs-corrected multiplier comparison. A from-scratch HPC run of the corrected Scripts 02–04 (with re-matching) is the exact version.

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

## Data archive and licence

- **Data archive:** the files in `data/` are not tracked in this repository (about 1.6 GB). They are deposited at Zenodo: **DOI: [TO BE ADDED ON DEPOSIT]**. Download and unpack the archive into `data/` before running any script.
- **Licence:** code is released under the MIT licence, see [LICENSE](LICENSE). Data and derived outputs in `data/` (Zenodo) are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); third-party inputs keep their original licences (see "Data not included" above).

## Computing environment

Figures were produced with R 4.6.1 on macOS (aarch64-apple-darwin23). Scripts 01–04 were run on an HPC cluster. The full `sessionInfo()` for the figure scripts is in [`sessionInfo.txt`](sessionInfo.txt). Package versions are listed in `sessionInfo.txt` but not pinned with `renv`. Scripts 05, 06 and 07 each run cleanly in a fresh R session (`Rscript scripts/<name>.R` from the package root); run them as separate processes, because sourcing them in one session lets attached packages mask each other (`tidyr::crossing`) and Script 06 then fails.

## Citation

If you use this code or data, please cite:

> Cho, F. H. T., Rhodes, J. R., Maron, M., Lunney, D., Iftekhar, M. S., McIlwee, A., Williams, B. A. (forthcoming). *Leakage effects of establishing protected areas in multi-use forests*

---

## Contact

For questions about the code or data, please contact the corresponding author or open an issue on the associated repository (frankie.cho@qut.edu.au)

---

## Contents of processed_data_for_plotting.RData

**Total objects:** 70 (data frames, lists, sf objects, ggplot objects, functions, scalars)

### Configuration / Paths

| Object | Type | Description |
|--------|------|-------------|
| `BASE_DIR` | `character` (length 1) | Base directory path for the project |
| `INPUT_DIR` | `character` (length 1) | Input data directory path |
| `OUTPUT_DIR` | `character` (length 1) | Output directory path |
| `scale` | `numeric` (length 1) | Numeric scale parameter |
| `ame_mdls_path` | `glue character` | Glue-formatted path to AME model outputs |

### Raw / Intermediate File Lists

These are data frames of file paths used to load model outputs, with columns `filename`, `run`, `year`.

| Object | Rows | Description |
|--------|------|-------------|
| `ame_all_matches` | 30 | File list for AME all-matches CSVs (buffer/neighbour × year) |
| `baci_model` | 29 | File list for BACI model outputs |
| `cum_clearing_all_matches` | 30 | File list for cumulative clearing all-matches CSVs |
| `match_data_pair` | 30 | File list for matched pair data |
| `match_data_pair_all` | 30 | File list for all matched pair data (across specs) |
| `matchit_summary` | 30 | File list for MatchIt summary outputs |
| `unmatch_data` | 30 | File list for unmatched data |

### Focal Forest (Treatment) Site Data

| Object | Class | Rows × Cols | Description |
|--------|-------|-------------|-------------|
| `focal_sf` | `sf` (MULTIPOLYGON) | 30 × 7 | Spatial polygons for focal (revoked) state forests. Columns: `run_id`, `woody_area_start_m2`, `woody_khab_area_start_m2`, `total_area_m2`, `woody_veg_lyr`, `five_year_clearing_ha`, `geom` |
| `focal_sf_tidy` | `data.frame` | 15 × 7 | Non-spatial tidy version of focal site attributes by year. Adds `woody_veg_area_ha` |
| `focal_sf_tidy_geom` | `sf` (MULTIPOLYGON) | 15 × 8 | Spatial version of `focal_sf_tidy` |
| `state_forests_revoked` | `sf` (MULTIPOLYGON) | 402 × 79 | Full attribute table for all revoked state forests. Columns include: `SFName`, `SFNo`, `revocation_year`, `area_km2`, `sf_id`, treated/untreated status, annual woody vegetation cover (`woody_veg_1988`–`woody_veg_2008`), annual SLATS clearing data (`slats_1988`–`slats_2020`), and covariates: `AgProf`, `DistCity`, `DistRoad`, `Drought`, `EcolCond`, `PopDen16`, `Eucalypt2018`, `Rainforest2018`, `Plantation2018`, `KoalaHabitatSuitability`, land use/tenure/zoning variables, `PolPref`, `PopGro16`, `remote2016`, plus geometry |

### Revocation Summaries

| Object | Class | Rows × Cols | Description |
|--------|-------|-------------|-------------|
| `revocation_summary` | `tbl_df` | 15 × 4 | Revocations per year: `revocation_year`, `n_revocations`, `revocation_area_ha`, `year` |
| `revocation_area_summary` | `tbl_df` | 30 × 5 | Area summaries by run and year: `run`, `year`, `n`, `woody_veg_focal_sf_ha`, `woody_veg_leakage_zone_ha` |
| `total_revocations` | `tbl_df` | 20 × 10 | Combined revocation + production data: `year`, `run`, `woody_area_start_m2`, `five_year_clearing_ha`, `total_area_ha`, `woody_veg_start_ha`, `revocation_year`, `n_revocations`, `revocation_area_ha`, `total_prod` |
| `sf_area_info` | `tbl_df` | 30 × 4 | Woody veg area summaries by run and year: `run`, `year`, `woody_veg_area_ha`, `woody_veg_leakage_zone_ha` |
| `woody_veg_area` | `tbl_df` | 30 × 7 | Woody vegetation area by run/year, linking to woody veg filenames and leakage zone area |

### Matching Results

| Object | Class | Rows × Cols | Description |
|--------|-------|-------------|-------------|
| `matchit_criteria` | `tbl_df` | 30 × 6 | Criteria flags per run/year: `ParallelTrends`, `SMD`, `PropKeep`, `meets_criteria` |
| `matchit_summary_df` | `tbl_df` | 573 × 14 | Summary statistics across matching specifications: method, distance metric, balance criteria flags (`ParallelTrends`, `SMD`, `PropKeep`), caliper, `model_selected`, `criteria_score`, `max_criteria_score`, `within_set`, `max_caliper` |
| `match_data_pair_df` | `tbl_df` | 30 × 12 | Aggregated BACI clearing estimates per run: `N`, `BC`, `AC`, `BI`, `AI`, `TOT_C`, `TOT_I`, `BACI`, `BACI_prop_C`, `BACI_prop_I` |
| `match_data_pair_all_summarised` | `tbl_df` | 136 × 13 | BACI estimates across all model specs (run_id × spec_id) |
| `match_data_pair_full` | `tbl_df` | 1,510,460 × 12 | Full plot-level matched panel data: `run_id`, `plot_id`, `treat`, `year`, `BA`, `CI`, `t`, `forest_clearing`, `prop_clearing`, `subclass`, `woody_veg_start` |
| `match_data_pair_subclass` | `tbl_df` | 55,287 × 13 | BACI estimates aggregated to matched-pair subclass level |
| `match_data_pair_summarised_range` | `tbl_df` | 16 × 10 | Median and range (glue-formatted) of BACI components by run/year |
| `match_data_pair_all_df` | `list` | — | List form of all matched pair data across specs |
| `sample_size` | `data.frame` | 29 × 2 | Sample size (`n`) per `run_id` |
| `woody_veg_ha_sample` | `data.frame` | 29 × 2 | Mean woody vegetation (ha) per `run_id` in matched sample |

### Cumulative Clearing Data

| Object | Class | Rows × Cols | Description |
|--------|-------|-------------|-------------|
| `match_data_cum_clearing` | `tbl_df` | 600 × 4 | Cumulative clearing time series for selected model: `run_id`, `CI` (confidence interval bound), `t` (time), `cum_clearing` |
| `cum_clearing_all_matches_df` | `tbl_df` | 24,000 × 19 | Cumulative clearing across all matching specs, including balance criteria flags and model selection indicators |

### Model Estimates (INLA / BACI)

| Object | Class | Rows × Cols / Length | Description |
|--------|-------|----------------------|-------------|
| `inla_mdls` | `list` | 29 elements | Fitted INLA model objects keyed by `run_year` (e.g. `buffer_1995`, `neighbour_2014`) |
| `baci_model_est` | `tbl_df` | 29 × 19 | INLA posterior estimates: `mean`, `sd`, quantiles (0.025, 0.05, 0.5, 0.95, 0.975, 0.01, 0.99), `mode`, `kld`, `sig_pos`, `sig_neg`, `N`, `N_treated`, `converge`, `all_nonzero` |
| `baci_model_ame_all_matches_df` | `grouped_df` | 29 × 18 | INLA estimates merged with site-level attributes and balance criteria |

### Average Marginal Effects (AME)

| Object | Class | Rows × Cols / Length | Description |
|--------|-------|----------------------|-------------|
| `ame_mdls` | `list` | 29 elements | AME model objects keyed by `run_year` |
| `ame_area` | `tbl_df` | 29 × 4 | AME estimates (2.5%, 50%, 97.5% quantiles) per `run_id` |
| `area_effect_df` | `grouped_df` | 29 × 32 | AME area effect with full site attributes, balance criteria, revocation metadata, and woody veg info |

### Leakage & TTE Estimates

| Object | Class | Rows × Cols | Description |
|--------|-------|-------------|-------------|
| `leakage_multiplier_df` | `grouped_df` | 19 × 37 | Primary leakage multiplier results. Includes AME quantiles, balance criteria, revocation info, and leakage/TTE-per-ha quantiles (`leak_*`, `tteph_*`) |
| `leakage_sf_tidy` | `sf` | 29 × 13 | Spatial leakage estimates by run/year: TTE (lb/median/ub), leakage (lb/median/ub), TTE-per-ha (lb/median/ub) |
| `pct_positive` | `tbl_df` | 2 × 2 | Proportion of simulations with positive leakage by run |
| `tte_dist` | `grouped_df` | 290,000 × 11 | Full posterior distribution of TTE estimates: `id`, `run`, `year`, `value`, `woody_veg_area_ha`, `woody_veg_leakage_zone_ha`, `tteph`, and balance criteria flags |
| `tte_year` | `grouped_df` | 29 × 5 | Cumulative TTE per hectare by run/year: `tteph_cum_lb`, `tteph_cum_median`, `tteph_cum_ub` |

### Distance / Grid Data

| Object | Class | Rows × Cols | Description |
|--------|-------|-------------|-------------|
| `avg_dist` | `grouped_df` | 20,000 × 3 | Average travel time to ephemeral harvest (`tteph`) per run and site `id` |
| `grid_points` | `tbl_df` | 712,714 × 2 | Grid of plot locations: `plot_id` and `geometry` (geoarrow format) |

### Production Data

| Object | Class | Rows × Cols | Description |
|--------|-------|-------------|-------------|
| `abares_prod` | `tbl_df` | 28 × 2 | ABARES annual timber production data: `Year` (1996–2023), `total_prod` |

### Saved ggplot Objects

| Object | Class | Description |
|--------|-------|-------------|
| `avg_dist_plt` | `ggplot` | Plot of average travel-time distances |
| `coef_plot` | `ggplot` | Coefficient/effect-size plot (likely BACI or AME estimates) |
| `cumulative_plt` | `ggplot` | Cumulative clearing over time |
| `cumulative_avg_dist_plt` | `patchwork ggplot` | Patchwork composite: cumulative clearing + distance |
| `total_revocations_plt` | `ggplot` | Plot of total revocations over time |

### Helper Functions

| Object | Description |
|--------|-------------|
| `add_unique_site_name` | Adds a unique site name column to a data frame |
| `data_non_zero` | Filters data to non-zero observations |
| `extract_slats_cov` | Extracts SLATS covariates from model data |
| `extract_years_from_files` | Parses years from file names |
| `generate_ame` | Generates average marginal effects |
| `generate_ame_irr` | Generates AME for IRR (incidence rate ratio) |
| `inla_is_converged` | Checks INLA model convergence |
| `plot_tag_theme` | ggplot theme helper for plot tags |
| `pull_ame_all_matches` | Loads AME all-matches CSVs from disk |
| `pull_cum_clearing` | Loads cumulative clearing CSVs |
| `pull_focal_sf` | Loads focal SF spatial data |
| `pull_inla_estimates` | Extracts INLA posterior estimates |
| `pull_inla_models` | Loads fitted INLA model objects |
| `pull_match_data_pair` | Loads matched pair data |
| `pull_match_data_pair_all` | Loads all-spec matched pair data |
| `pull_matchit_summary` | Loads MatchIt balance summaries |
| `rename_run` | Renames run identifiers for display |
| `summarise_cum_clearing` | Summarises cumulative clearing data |
| `summarise_match_data` | Summarises match-level data |
| `summarise_match_data_plot_level` | Summarises match data at plot level |
| `treatment_control_plot` | Generates treatment vs. control diagnostic plots |

### Key Run IDs

Models are estimated separately for two matching strategies × 15 cohort years:

- **`buffer_YYYY`** — buffer-zone matching (YYYY ∈ {1995–2014})
- **`neighbour_YYYY`** — neighbour matching (YYYY ∈ {1995–2014})

29 combinations total (one neighbour year differs: includes 2012 but not 2013).
