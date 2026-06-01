# =============================================================================
# Script 01: Extract pixel-level covariates and deforestation data
# =============================================================================
# PURPOSE
#   Extracts annual deforestation (SLATS), woody vegetation cover, protected-
#   area coverage, and spatial covariates for each 10 km grid cell covering
#   New South Wales (NSW). Outputs a single Parquet file used in Script 02.
#
# REFERENCE SCRIPT
#   This script was run on an HPC cluster with large raster inputs not
#   included in this package (SLATS tiles, covariate rasters, woody veg
#   layers, ~100 GB total). It is provided for methodological transparency.
#   To replicate figures, start from Script 05 with the pre-computed
#   data/processed_data_for_plotting.RData.
#
# INPUTS (HPC paths)
#   {INPUT_DIR}/resampled/slats_resampled.tif
#   {INPUT_DIR}/woody_veg/annual/woody_veg_*.tif          (annual layers)
#   {INPUT_DIR}/protected_area_years/PA_*.tif             (annual layers)
#   {INPUT_DIR}/covariates_raster/*.tif                   (25+ rasters)
#   {INPUT_DIR}/nsw_grid/nsw_grid_points.shp              (generated here)
#
# OUTPUT
#   {INPUT_DIR}/matching_data/parquet/matching_data_pixel_{job_id}.parquet
#   (or matching_data_pixel.parquet when run locally for a small sample)
# =============================================================================

library(dplyr)
library(sf)
library(exactextractr)
library(arrow)
library(glue)
library(terra)

# HPC job array index (set by PBS scheduler; absent for local runs)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && !is.na(args[1])) job_id <- as.numeric(args[1])

if (exists("job_id")) {
  home      <- "/mnt/hpccs01/home/chof"
} else {
  home      <- "H:"
}

INPUT_DIR  <- file.path(home, "forestry-leakage/input")
OUTPUT_DIR <- file.path(home, "forestry-leakage/output")

use_woody_veg_annual <- TRUE

# -- SLATS deforestation rasters (resampled to 25 m) --------------------------
slats_name <- list.files(
  file.path(home, "nsw-habitat-loss/output/slats"),
  pattern = "*.tif$", full.names = FALSE
) %>%
  sapply(basename) %>%
  sapply(tools::file_path_sans_ext) %>%
  unname()

slats_rast_resampled <- rast(file.path(INPUT_DIR, "resampled/slats_resampled.tif"))
names(slats_rast_resampled) <- slats_name

# -- Woody vegetation layers --------------------------------------------------
woody_veg <- rast(file.path(INPUT_DIR, "woody_veg/woody_nsw.tif"))

if (use_woody_veg_annual) {
  woody_veg_lyrs <- list.files(
    file.path(INPUT_DIR, "woody_veg/annual"),
    pattern = "woody_veg.*\\.tif$", full.names = TRUE
  )
  woody_veg_ann  <- rast(woody_veg_lyrs[1:13])
  names(woody_veg_ann) <- tools::file_path_sans_ext(basename(woody_veg_lyrs[1:13]))
}

# -- Protected area layers (one per year) -------------------------------------
pa_lyrs     <- list.files(file.path(INPUT_DIR, "protected_area_years"), pattern = "PA_.*\\.tif$", full.names = TRUE)
pa_lyrs_ann <- rast(pa_lyrs)
names(pa_lyrs_ann) <- tools::file_path_sans_ext(basename(pa_lyrs))

# -- Covariate rasters --------------------------------------------------------
# Helper to expand multi-layer rasters (e.g. PCA outputs) into named layers
insert_items <- function(vec, items, pos) {
  c(vec[seq_len(pos - 1)], items, vec[(pos + 1):length(vec)])
}

covariate_list   <- list.files(file.path(INPUT_DIR, "covariates_raster"), pattern = "*.tif$", full.names = TRUE)
covariates       <- rast(covariate_list)
covariate_names  <- tools::file_path_sans_ext(basename(covariate_list))
covariate_names  <- insert_items(covariate_names, c("ScEc_PC1","ScEc_PC2","ScEc_PC3","ScEc_PC4","ScEc_PC5"),
                                 which(covariate_names == "SocioEcon16_PC"))
covariate_names  <- insert_items(covariate_names, c("Soil_PC1","Soil_PC2","Soil_PC3"),
                                 which(covariate_names == "TSoilPC"))
names(covariates) <- covariate_names

cont_covs <- c("AgProf","DistCity","DistRoad","Drought","EcolCond","PopDen16",
               "Eucalypt2018","Rainforest2018","Plantation2018",
               "temp","prec","slope",
               "ScEc_PC1","ScEc_PC2","ScEc_PC3","ScEc_PC4","ScEc_PC5",
               "Soil_PC1","Soil_PC2","Soil_PC3")
disc_covs <- c("KoalaHabitatSuitability","LandTen","LandUse","LandZone","NatVegReg",
               "NSW_forten18_ForCode","NSW_forten18_ForTen","NSW_forten18_ForTenType",
               "NSW_forten18_ForType","NSW_forten18_TenType",
               "PlanZone","PolPref","Fire","PopGro16","remote2016")

# -- NSW grid (10 km regular grid; generated on first run, cached thereafter) -
nsw_border         <- st_read(file.path(INPUT_DIR, "nsw_border/NSW_LGA_GDA2020/NSW_LGA_GDA2020.shp")) %>% st_union()
grid_points_path   <- file.path(INPUT_DIR, "nsw_grid/nsw_grid_points.shp")

if (file.exists(grid_points_path)) {
  grid_points_proj <- st_read(grid_points_path)
} else {
  grid            <- st_make_grid(nsw_border, cellsize = 10000)
  grid_points     <- st_centroid(grid) %>% st_sfc() %>% st_sf()
  grid_points     <- grid_points[lengths(st_intersects(grid_points, nsw_border)) > 0, ]
  grid_points_proj <- grid_points %>%
    st_buffer(150) %>%
    vect() %>%
    project(covariates) %>%
    st_as_sf() %>%
    mutate(plot_id = row_number())
  grid_points_proj <- grid_points_proj[!st_is_empty(grid_points_proj), ]
  st_write(grid_points_proj, grid_points_path, append = FALSE)
}

# -- Subset grid for HPC job array --------------------------------------------
if (exists("job_id")) {
  n_iterations   <- 1000
  set.seed(348945)
  n_points       <- nrow(grid_points_proj)
  chunk_size     <- ceiling(n_points / n_iterations)
  selection_seq  <- (job_id - 1) * chunk_size + seq_len(chunk_size)
  selection_seq  <- selection_seq[selection_seq <= n_points]
  grid_points_sample <- grid_points_proj[selection_seq, ]
} else {
  grid_points_sample <- grid_points_proj
}

# -- Extract values at each grid point ----------------------------------------
if (use_woody_veg_annual) {
  woody_extract <- exact_extract(woody_veg_ann, grid_points_sample, coverage_area = TRUE,
    fun = function(x, coverage_area) {
      bind_cols(x, coverage_area = coverage_area) %>%
        mutate_at(names(woody_veg_ann), ~ if_else(is.na(.x) | (.x < 1), 0, (.x >= 1) * coverage_area)) %>%
        colSums(na.rm = TRUE)
    }, force_df = TRUE)

  woody_values <- woody_extract %>% t() %>% as.data.frame()
  rownames(woody_values) <- NULL
  colnames(woody_values)  <- names(woody_veg_ann)
}

slats_extract <- exact_extract(slats_rast_resampled, grid_points_sample, coverage_area = TRUE,
  fun = function(x, coverage_area) {
    bind_cols(x, coverage_area = coverage_area) %>%
      mutate_at(names(slats_rast_resampled), ~ if_else(is.na(.x) | (.x != 4), 0, (.x == 4) * coverage_area)) %>%
      colSums(na.rm = TRUE)
  }, force_df = TRUE)

slats_values <- slats_extract %>% t() %>% as.data.frame() %>% dplyr::select(-coverage_area)
rownames(slats_values) <- NULL

# SLATS data contains bi-annual observations pre-2006; split evenly across years
biannual_deforestation        <- slats_values[, rep(1:9, each = 2)] / 2
names(biannual_deforestation) <- paste0("slats_", 1988:2005)
annual_deforestation          <- slats_values[, 10:ncol(slats_values)]
names(annual_deforestation)   <- paste0("slats_", 2006:2020)
ann_deforestation_area        <- bind_cols(biannual_deforestation, annual_deforestation)

cont_cov_extract <- exact_extract(covariates[[cont_covs]], grid_points_sample, fun = "mean") %>%
  rename_all(~ stringr::str_replace_all(.x, "mean.", ""))
disc_cov_extract <- exact_extract(covariates[[disc_covs]], grid_points_sample, fun = "mode") %>%
  rename_all(~ stringr::str_replace_all(.x, "mode.", ""))

pa_values <- exact_extract(pa_lyrs_ann, grid_points_sample, fun = "mean", default_value = 0) %>%
  rename_all(~ stringr::str_replace_all(.x, "mean.", ""))

# -- Bind and write -----------------------------------------------------------
if (use_woody_veg_annual) {
  df <- bind_cols(plot_id = grid_points_sample$plot_id,
                  woody_values, ann_deforestation_area,
                  cont_cov_extract, disc_cov_extract, pa_values)
} else {
  df <- bind_cols(plot_id = grid_points_sample$plot_id,
                  ann_deforestation_area, cont_cov_extract, disc_cov_extract, pa_values)
}

if (exists("job_id")) {
  write_parquet(df, file.path(INPUT_DIR, glue("matching_data/parquet/matching_data_pixel_{job_id}.parquet")))
} else {
  write_parquet(df, file.path(INPUT_DIR, "matching_data/matching_data_pixel.parquet"))
}
