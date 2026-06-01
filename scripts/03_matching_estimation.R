# =============================================================================
# Script 03: BACI estimation via INLA negative-binomial models
# =============================================================================
# PURPOSE
#   For each revocation year (and leakage channel), fits a Bayesian Before-
#   After Control-Impact (BACI) model using INLA with a negative-binomial
#   family. The BACI coefficient captures the differential change in forest
#   clearing between treated and matched control plots.
#
#   Also generates Average Marginal Effect (AME) estimates in hectares, and
#   runs a distance-threshold sensitivity analysis (1–10 km from the focal SF).
#
# REFERENCE SCRIPT
#   Requires model-selection outputs from Script 02. Run on HPC.
#
# INPUTS (HPC paths)
#   {OUTPUT_DIR}/match_data_pair_{buffer,neighbour}/match_data_long_*.parquet
#   {OUTPUT_DIR}/matchit_{buffer,neighbour}/matchit_*_summary_*.csv
#   {OUTPUT_DIR}/focal_state_forest_gpkg/focal_sf_year.gpkg
#   {INPUT_DIR}/nsw_grid/nsw_grid_points.parquet
#
# OUTPUTS
#   {OUTPUT_DIR}/baci_models_{buffer,neighbour}/baci_model_*_{site_group}.rds
#   {OUTPUT_DIR}/ame/baci_model_ame_*_{site_group}.csv
#   {OUTPUT_DIR}/ame/dist/baci_model_ame_*_{site_group}_{d}.csv
#   {OUTPUT_DIR}/plot_data/matching_df_{run}.csv
#   {OUTPUT_DIR}/plot_data/matchit_simple_summary_{run}.csv
# =============================================================================

library(tidyverse)
library(future)
library(future.apply)
library(glue)
library(arrow)
library(MatchIt)
library(INLA)
library(sf)
library(geoarrow)

source("scripts/ame.R")

INPUT_DIR  <- "H:/forestry-leakage/input"
OUTPUT_DIR <- "H:/forestry-leakage/output/run_202508"

grid_points <- open_dataset(file.path(INPUT_DIR, "nsw_grid/nsw_grid_points.parquet"))

run_site           <- FALSE   # Estimate individual site-level effects
run_dist_threshold <- TRUE    # Run 1–10 km distance sensitivity analysis

# ── Helper: extract and sort file names by numeric site ID ───────────────────
sort_by_id <- function(x, i = -2) {
  id <- as.numeric(sapply(tools::file_path_sans_ext(x), str_split_i, pattern = "_", i = i))
  tools::file_path_sans_ext(x[sort.int(id, index.return = TRUE)$ix])
}

get_id <- function(x, i = -2) {
  as.numeric(sapply(tools::file_path_sans_ext(x), str_split_i, pattern = "_", i = i))
}

filter_by_id <- function(x, filter_id = 72:86, ...) {
  x[as.character(get_id(x, ...)) %in% as.character(filter_id)]
}

extract_site_id_from_filename <- function(x) sapply(x, str_split_i, pattern = "_", i = -2) %>% unname()
extract_name_from_filename    <- function(x) sapply(x, str_split_i, pattern = "_", i = -1) %>% unname()

# ── Main estimation loop (buffer and neighbour channels) ─────────────────────
for (run in c("buffer", "neighbour")) {

  match_data_pair_path       <- glue("{OUTPUT_DIR}/match_data_pair_{run}")
  matchit_path               <- glue("{OUTPUT_DIR}/matchit_{run}")
  unmatch_data_path          <- glue("{OUTPUT_DIR}/unmatch_data_{run}")
  treated_gpkg_path          <- glue("{OUTPUT_DIR}/treated_polygon_gpkg")
  focal_state_forest_gpkg_path <- glue("{OUTPUT_DIR}/focal_state_forest_gpkg")

  match_data_pair_filename <- list.files(match_data_pair_path, pattern = "\\.parquet$") %>% sort_by_id() %>% filter_by_id()
  match_data_all_filename  <- list.files(match_data_pair_path, pattern = "\\.rds$")     %>% sort_by_id() %>% filter_by_id()
  match_summary_filename   <- list.files(matchit_path,         pattern = "\\.csv$")     %>% sort_by_id() %>% filter_by_id()
  focal_sf_gpkg_filename   <- list.files(focal_state_forest_gpkg_path, pattern = glue("^focal_state_forest_{run}")) %>%
    sort_by_id() %>% filter_by_id()

  focal_sf_gpkg_year <- st_read(glue("{OUTPUT_DIR}/focal_state_forest_gpkg/focal_sf_year.gpkg"))

  site_names    <- extract_name_from_filename(match_data_pair_filename)
  site_id       <- extract_site_id_from_filename(match_data_pair_filename)
  site_name_fmt <- as.list(glue("{site_id}. {str_to_title(site_names)}")) %>% unlist()

  # Load MatchIt summary tables
  matchit_summary_list <- lapply(match_summary_filename, function(name) {
    read_csv(glue("{matchit_path}/{tools::file_path_sans_ext(name)}.csv"), show_col_types = FALSE)
  })
  names(matchit_summary_list) <- site_name_fmt

  matchit_summary_df <- matchit_summary_list %>%
    lapply(\(x) { x$covariate <- as.character(x$covariate); x }) %>%
    bind_rows(.id = "site_group")

  # Load matched pair data
  match_data_pair_list <- lapply(match_data_pair_filename, function(name) {
    read_parquet(glue("{match_data_pair_path}/{tools::file_path_sans_ext(name)}.parquet"),
                 col_select = c("plot_id","treat","year","BA","CI","t",
                                "forest_clearing","prop_clearing","subclass",
                                "woody_veg_start","KoalaHabitatSuitability"))
  })
  names(match_data_pair_list) <- site_name_fmt

  matching_df <- match_data_pair_list %>%
    bind_rows(.id = "site_group") %>%
    mutate(site    = str_split_i(site_group, ". ", 2),
           site_id = str_split_i(site_group, ". ", 1)) %>%
    arrange(as.numeric(site_id))

  site_names_df <- data.frame(site_id = site_id, site_name = site_names, site_group = site_name_fmt)
  site_names_df$site_group <- factor(site_names_df$site_group, site_names_df$site_group)

  # Save matching data and summary
  write_csv(matching_df,                           glue("{OUTPUT_DIR}/plot_data/matching_df_{run}.csv"))
  write_csv(
    matchit_summary_df %>%
      filter(model_selected) %>%
      group_by(site_group, distance, link, replace, caliper, new_vars) %>%
      mutate(meets_balance_criteria = all(ParallelTrends, SMD, PropKeep)) %>%
      summarise(spec_id = first(spec_id), meets_balance_criteria = any(meets_balance_criteria),
                parallel_trends = any(ParallelTrends), SMD = any(SMD), PropKeep = any(PropKeep)),
    glue("{OUTPUT_DIR}/plot_data/matchit_simple_summary_{run}.csv")
  )

  baci_model_dir <- glue("{OUTPUT_DIR}/baci_models_{run}")
  if (!dir.exists(baci_model_dir)) dir.create(baci_model_dir)

  plan(sequential)

  for (i in seq_along(unique(matching_df$site_group))) {
    site_group <- unique(matching_df$site_group)[i]

    df <- matching_df[matching_df$site_group == site_group, ] %>%
      mutate(forest_clearing_cell = ceiling(forest_clearing / 625),
             total_woody_cells    = ceiling(woody_veg_start / 625),
             BACI = BA * CI) %>%
      group_by(plot_id) %>% arrange(plot_id, t) %>% ungroup()

    if (sum(df$forest_clearing_cell) == 0) {
      message(site_group, ": skipped — no clearing in observation period")
      next()
    }
    if (!run_site && !grepl("\\d{4,}", site_group)) {
      message(site_group, ": skipped — run_site = FALSE")
      next()
    }

    df_summarised <- df %>%
      group_by(site_group, plot_id, BA, CI, subclass) %>%
      summarise(forest_clearing_cell = sum(forest_clearing_cell),
                total_woody_cells    = max(total_woody_cells), .groups = "drop") %>%
      rowwise() %>%
      mutate(forest_clearing_cell = min(forest_clearing_cell, total_woody_cells),
             BACI = BA * CI)

    baci_formula <- forest_clearing_cell ~ 1 + BA + CI + BACI + offset(log(total_woody_cells))

    v <- inla(baci_formula, data = df_summarised, family = "nbinomial",
              control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE),
              verbose = FALSE)
    write_rds(v, file.path(baci_model_dir, glue("baci_model_{run}_{site_group}.rds")), compress = "none")

    ame <- generate_ame_irr(v, n_samples = 10000)
    write_csv(as.data.frame(ame), glue("{OUTPUT_DIR}/ame/baci_model_ame_{run}_{site_group}.csv"))

    # Distance-threshold sensitivity: re-estimate restricting treatment to plots ≤ d km from focal SF
    if (grepl("\\d{4,}", site_group) && run_dist_threshold) {
      year    <- substr(site_group, nchar(site_group) - 3, nchar(site_group))
      focal_sf <- focal_sf_gpkg_year[
        year == substr(focal_sf_gpkg_year$run_id, nchar(focal_sf_gpkg_year$run_id) - 3, nchar(focal_sf_gpkg_year$run_id)),
      ][1, ]

      grid_points_loaded <- grid_points %>%
        filter(plot_id %in% df_summarised$plot_id) %>% st_as_sf() %>% st_centroid()

      dist_to_sf <- grid_points_loaded %>%
        filter(plot_id %in% df_summarised$plot_id[df_summarised$CI == 1]) %>%
        mutate(dist = as.numeric(units::set_units(st_distance(., focal_sf), "km"))) %>%
        st_drop_geometry()

      dist_dir <- file.path(baci_model_dir, "dist")
      if (!dir.exists(dist_dir)) dir.create(dist_dir)

      for (d in 1:10) {
        within_dist_plotid <- dist_to_sf %>% filter(dist < d) %>% pull(plot_id)
        dist_subclasses    <- df_summarised$subclass[df_summarised$plot_id %in% within_dist_plotid]
        df_dist            <- filter(df_summarised, subclass %in% dist_subclasses)

        v_dist <- inla(baci_formula, data = df_dist, family = "nbinomial",
                       control.compute = list(dic = TRUE, cpo = TRUE, waic = TRUE, config = TRUE),
                       verbose = FALSE)
        write_rds(v_dist, file.path(dist_dir, glue("baci_model_{run}_{site_group}_{d}.rds")), compress = "none")

        ame_dist <- generate_ame_irr(v_dist, n_samples = 10000)
        write_csv(as.data.frame(ame_dist), glue("{OUTPUT_DIR}/ame/dist/baci_model_ame_{run}_{site_group}_{d}.csv"))
      }
    }
  }

  plan(sequential)
}
