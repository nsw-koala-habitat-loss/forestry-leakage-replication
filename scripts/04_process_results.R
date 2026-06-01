# =============================================================================
# Script 04: Aggregate model outputs into a single processed data file
# =============================================================================
# PURPOSE
#   Loads all INLA model objects, MatchIt summaries, and matched pair data
#   produced by Scripts 02–03. Computes:
#     • BACI coefficient estimates and credible intervals
#     • Average Marginal Effects (AME) in hectares
#     • Leakage multipliers and treatment effect per hectare of protected area
#     • Revocation and log-production summaries for plotting
#
#   Saves everything needed by Script 05 into a single .RData file.
#
# REFERENCE SCRIPT
#   Requires INLA model RDS files from Script 03. Run locally on the PC with
#   full output data. The pre-computed output (data/processed_data_for_plotting.RData)
#   is included in this package so you can skip directly to Script 05.
#
# INPUTS (relative to project root or HPC paths)
#   {OUTPUT_DIR}/baci_models_{buffer,neighbour}/*.rds
#   {OUTPUT_DIR}/matchit_{buffer,neighbour}/*.csv
#   {OUTPUT_DIR}/match_data_pair_{buffer,neighbour}/*.parquet
#   {OUTPUT_DIR}/match_data_pair_{buffer,neighbour}/*.rds  (all-specs)
#   {OUTPUT_DIR}/cum_clearing_all_matches/*.csv
#   {OUTPUT_DIR}/ame/*.csv
#   {OUTPUT_DIR}/focal_state_forest_gpkg/focal_sf_year.gpkg
#   {OUTPUT_DIR}/plot_data/year_run_df.csv
#   {INPUT_DIR}/state_forests/state_forests_revoked_update.gpkg
#   {INPUT_DIR}/nsw_grid/nsw_grid_points.parquet
#   ext_data/*.xlsx                                         (ABARES log production)
#
# OUTPUT
#   data/processed_data_for_plotting.RData
# =============================================================================

library(ggplot2)
library(tidyverse)
library(glue)
library(sf)
library(purrr)
library(readr)
library(INLA)
library(arrow)
library(geoarrow)

source("scripts/config.R")
source("scripts/helpers.R")
source("scripts/ame.R")

# Paths to the full model output (adjust to your local/HPC path)
OUTPUT_DIR_MODELS <- "H:/forestry-leakage/output/run_202508"
INPUT_DIR <- "H:/forestry-leakage/input"

grid_points <- read_parquet(file.path(
  INPUT_DIR,
  "nsw_grid/nsw_grid_points.parquet"
))

# ── Enumerate model output files ─────────────────────────────────────────────

# Default pattern matches:  baci_model_{run}_{nn}. Y{year}.rds
baci_model <- c(
  list.files(glue("{OUTPUT_DIR_MODELS}/baci_models_buffer"), full.names = TRUE),
  list.files(
    glue("{OUTPUT_DIR_MODELS}/baci_models_neighbour"),
    full.names = TRUE
  )
) %>%
  extract_years_from_files()

# Override pattern matches:  unmatch_data_{run}_{nn}_Y{year}.rds  (no space before Y)
unmatch_data <- c(
  list.files(
    glue("{OUTPUT_DIR_MODELS}/unmatch_data_buffer"),
    full.names = TRUE
  ),
  list.files(
    glue("{OUTPUT_DIR_MODELS}/unmatch_data_neighbour"),
    full.names = TRUE
  )
) %>%
  extract_years_from_files(pattern = "_([[:alnum:]]+)_\\d+\\_Y(\\d{4})")

woody_veg_area <- read_csv(glue(
  "{OUTPUT_DIR_MODELS}/plot_data/year_run_df.csv"
)) %>%
  mutate(woody_veg_leakage_zone_ha = woody_veg_area * 0.0001)

# Pattern matches the selected-spec files:  match_data_long_{nn}_Y{year}.parquet
match_data_pair <- c(
  list.files(
    glue("{OUTPUT_DIR_MODELS}/match_data_pair_buffer"),
    full.names = TRUE
  ),
  list.files(
    glue("{OUTPUT_DIR_MODELS}/match_data_pair_neighbour"),
    full.names = TRUE
  )
) %>%
  extract_years_from_files(
    pattern = "_([[:alnum:]]+)/match_data_long_\\d+\\_Y(\\d{4})",
    extension = ".parquet"
  )

# Pattern matches the all-specs files:  match_data_long_all_{nn}_Y{year}.rds
match_data_pair_all <- c(
  list.files(
    glue("{OUTPUT_DIR_MODELS}/match_data_pair_buffer"),
    full.names = TRUE
  ),
  list.files(
    glue("{OUTPUT_DIR_MODELS}/match_data_pair_neighbour"),
    full.names = TRUE
  )
) %>%
  extract_years_from_files(
    pattern = "_([[:alnum:]]+)/match_data_long_all_\\d+\\_Y(\\d{4})",
    extension = ".rds"
  )

matchit_summary <- c(
  list.files(glue("{OUTPUT_DIR_MODELS}/matchit_buffer"), full.names = TRUE),
  list.files(glue("{OUTPUT_DIR_MODELS}/matchit_neighbour"), full.names = TRUE)
) %>%
  extract_years_from_files(
    pattern = "_([[:alnum:]]+)_summary_\\d+\\_Y(\\d{4})",
    extension = ".csv"
  )

cum_clearing_all_matches <- list.files(
  glue("{OUTPUT_DIR_MODELS}/cum_clearing_all_matches"),
  full.names = TRUE
) %>%
  extract_years_from_files(
    pattern = "_([[:alnum:]]+)_\\d+\\_Y(\\d{4})",
    extension = ".csv"
  )

# ── Load revocation shapefile and ABARES log-production data ──────────────────

if (!exists("state_forests_revoked")) {
  state_forests_revoked <- st_read(
    file.path(DATA_DIR, "state_forests/state_forests_revoked_update.gpkg"),
    quiet = TRUE
  )
}

abares_xlsx <- file.path(
  "ext_data",
  "5 Dashboard data tables - Plantation area and log production AFWPS Production and Trade to 2023-24 – MS Excel.xlsx"
)
if (file.exists(abares_xlsx)) {
  abares_prod <- readxl::read_excel(
    abares_xlsx,
    sheet = "Volume of production"
  ) %>%
    filter(State == "New South Wales" & Category_1 == "Log production") %>%
    mutate(m3 = if_else(Units == "’000  m3", Value * 1000, Value)) %>%
    group_by(Year, Category_2) %>%
    summarise(total_prod = sum(m3))
} else {
  warning(
    "ABARES Excel file not found; log production data will be absent from plots."
  )
  abares_prod <- tibble(
    Year = integer(),
    Category_2 = character(),
    total_prod = numeric()
  )
}

# ── Revocation summary ────────────────────────────────────────────────────────

revocation_summary <- state_forests_revoked %>%
  mutate(area_ha = as.numeric(units::set_units(st_area(.), "ha"))) %>%
  st_drop_geometry() %>%
  group_by(revocation_year) %>%
  summarise(n_revocations = n(), revocation_area_ha = sum(area_ha)) %>%
  mutate(year = revocation_year)

# ── Load INLA models and compute AME ─────────────────────────────────────────

inla_mdls <- pull_inla_models(baci_model)
names(inla_mdls) <- paste0(baci_model$run, "_", baci_model$year)

ame_mdls_path <- file.path(OUTPUT_DIR_MODELS, "ame/ame_mdls.rds")
if (!file.exists(ame_mdls_path)) {
  ame_mdls <- map(
    inla_mdls,
    \(x) generate_ame_irr(x, n_samples = 10000),
    .progress = "Drawing AME"
  )
  write_rds(ame_mdls, ame_mdls_path)
} else {
  ame_mdls <- read_rds(ame_mdls_path)
}

ame_area <- ame_mdls %>%
  map(pull, var = "delta_area") %>%
  bind_rows() %>%
  apply(2, quantile, p = c(.025, .5, .975), na.rm = TRUE) %>%
  t() %>%
  as_tibble() %>%
  mutate(across(everything(), \(x) round(x, 2)), run_id = names(inla_mdls))

woody_veg_ha_sample <- map_dbl(inla_mdls, \(x) {
  sum(x$.args$data$total_woody_cells) * 625 * 0.0001
}) %>%
  as.data.frame() %>%
  rownames_to_column("run_id") %>%
  rename(woody_veg_ha = ".")

# ── MatchIt summary processing ────────────────────────────────────────────────

matchit_summary_df <- pull_matchit_summary(matchit_summary) %>%
  group_by(run_id, spec_id) %>%
  summarise(across(
    c(
      method,
      distance,
      ParallelTrends,
      SMD,
      PropKeep,
      replace,
      caliper,
      model_selected
    ),
    first
  )) %>%
  mutate(criteria_score = ParallelTrends + SMD + PropKeep) %>%
  group_by(run_id) %>%
  mutate(max_criteria_score = max(criteria_score, na.rm = TRUE)) %>%
  mutate(
    within_set = (criteria_score == max_criteria_score) | model_selected
  ) %>%
  mutate(max_caliper = max(caliper, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    within_set = if_else(
      is.na(caliper) | caliper == max_caliper | model_selected,
      within_set,
      FALSE
    )
  )

matchit_criteria <- matchit_summary_df %>%
  group_by(run_id) %>%
  filter(model_selected) %>%
  summarise(
    ParallelTrends = all(ParallelTrends),
    SMD = all(SMD),
    PropKeep = all(PropKeep)
  ) %>%
  separate(run_id, c("run", "year"), sep = "_") %>%
  mutate(meets_criteria = (ParallelTrends + SMD + PropKeep) == 3)

baci_model_est <- pull_inla_estimates(baci_model) %>%
  separate(run_id, c("run", "year"), sep = "_")

sample_size <- map_dbl(inla_mdls, \(x) {
  length(unique(x$.args$data$plot_id[x$.args$data$CI == 1]))
}) %>%
  as.data.frame() %>%
  rownames_to_column("run_id") %>%
  rename(n = ".")

# ── Focal state forest shapefile ──────────────────────────────────────────────

focal_sf <- st_read(file.path(
  OUTPUT_DIR_MODELS,
  "focal_state_forest_gpkg/focal_sf_year.gpkg"
))
focal_sf_tidy <- st_drop_geometry(focal_sf) %>%
  mutate(woody_veg_area_ha = woody_area_start_m2 * 0.0001) %>%
  distinct(run_id, .keep_all = TRUE) %>%
  separate(run_id, c("run", "year")) %>%
  dplyr::select(-run)

# ── Proportional AME and area effect ─────────────────────────────────────────

# Note: `2.5%`/`50%`/`97.5%` are the raw AME credible interval bounds in hectares.
# `0.025quant`/`0.5quant`/`0.975quant` are the same values divided by woody_veg_area_ha
# (TTE per ha of protected area). Both naming conventions are carried forward.
baci_model_ame_all_matches_df <- ame_area %>%
  mutate(id = row_number()) %>%
  left_join(woody_veg_ha_sample) %>%
  separate(run_id, c("run", "year"), sep = "_") %>%
  left_join(focal_sf_tidy) %>%
  group_by(run, year) %>%
  mutate(
    `0.025quant` = `2.5%` / woody_veg_area_ha,
    `0.5quant` = `50%` / woody_veg_area_ha,
    `0.975quant` = `97.5%` / woody_veg_area_ha
  ) %>%
  mutate(
    sig_pos = `0.975quant` > 0 & `0.025quant` > 0 & `0.5quant` > 0,
    sig_neg = `0.975quant` < 0 & `0.025quant` < 0 & `0.5quant` < 0
  )

match_data_pair_full <- pull_match_data_pair(match_data_pair)
match_data_pair_subclass <- pull_match_data_pair(
  match_data_pair,
  sum_func = summarise_match_data_plot_level
)
match_data_pair_df <- pull_match_data_pair(
  match_data_pair,
  sum_func = summarise_match_data
)
match_data_cum_clearing <- pull_match_data_pair(
  match_data_pair,
  sum_func = summarise_cum_clearing
)

cum_clearing_all_matches_df <- pull_cum_clearing(cum_clearing_all_matches) %>%
  left_join(matchit_summary_df, by = c("run_id", "spec_id"))

match_data_pair_all_df <- pull_match_data_pair_all(
  match_data_pair_all,
  matchit_summary_df
)
match_data_pair_all_summarised <- match_data_pair_all_df %>%
  bind_rows(.id = "run_id") %>%
  summarise_match_data(summary_vars = c("run_id", "spec_id"))

# ── Area effect and leakage multiplier ───────────────────────────────────────

area_effect_df <- baci_model_ame_all_matches_df %>%
  left_join(matchit_criteria) %>%
  left_join(dplyr::select(
    baci_model_est,
    c("run", "year", "all_nonzero", "N_treated")
  )) %>%
  mutate(year = as.numeric(year)) %>%
  left_join(revocation_summary) %>%
  left_join(woody_veg_area, by = c("run", "year")) %>%
  mutate(
    `0.5quant` = `50%` * woody_veg_leakage_zone_ha / woody_veg_area_ha,
    `0.025quant` = `2.5%` * woody_veg_leakage_zone_ha / woody_veg_area_ha,
    `0.975quant` = `97.5%` * woody_veg_leakage_zone_ha / woody_veg_area_ha
  ) %>%
  mutate(year = as.character(year))

# leakage_multiplier = TTE_ha / five_year_clearing_ha, where five_year_clearing_ha
# is observed clearing *inside* the focal SF in the 5 years before revocation.
# A value of 1 (100%) means leakage fully offsets the conservation gain.
# tteph = TTE_ha / woody_veg_area_ha: treatment effect per hectare of new protected area.
leakage_multiplier_df <- area_effect_df %>%
  dplyr::select(-any_of("woody_veg_area")) %>%
  left_join(focal_sf_tidy) %>%
  mutate(
    `leak_0.025quant` = `0.025quant` / five_year_clearing_ha,
    `leak_0.5quant` = `0.5quant` / five_year_clearing_ha,
    `leak_0.975quant` = `0.975quant` / five_year_clearing_ha,
    `tteph_0.025quant` = `0.025quant` / woody_veg_area_ha,
    `tteph_0.5quant` = `0.5quant` / woody_veg_area_ha,
    `tteph_0.975quant` = `0.975quant` / woody_veg_area_ha
  ) %>%
  rename_run() %>%
  filter(all_nonzero) %>%
  mutate(year = as.numeric(year))

# ── Save all objects for Script 05 ───────────────────────────────────────────

save(
  list = ls(),
  file = file.path(DATA_DIR, "processed_data_for_plotting.RData")
)
message(
  "Saved processed data to: ",
  file.path(DATA_DIR, "processed_data_for_plotting.RData")
)
