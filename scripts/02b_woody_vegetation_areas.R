# =============================================================================
# Script 02b: Woody-vegetation areas and on-site avoided clearing
# =============================================================================
# Verbatim port of scripts/woody_veg_calc.R. Per focal SF (buffer channel),
# from the annual woody-veg + SLATS rasters, computes:
#   woody_area_start_m2    woody veg in the focal SF at the revocation year
#   five_year_clearing_ha  on-site avoided clearing = SLATS class-4 clearing,
#                          sum(pre[y-5..y-1]) - sum(post[y+1..y+5]), y excluded
#   woody_veg_area         woody veg in the 10 km treated belt -> year_run_df.csv
#
# Runs between Scripts 02 and 03. Reference script: needs the upstream rasters;
# focal_sf_year.gpkg is shipped pre-computed in data/.
#
# Human-source bugs, preserved (not fixed):
#   (1) `unmatch_data` lists unmatch_data_buffer twice, neighbour never -> the
#       gpkg holds only buffer_{year} (each duplicated). Benign: Scripts 04/05
#       distinct(run_id) and join by year, so neighbour reuses the same focal value.
#   (2) ann_deforestation_area, joined_df$year.x and focal_sf are aligned by row
#       position across independently-sorted lists, not by a key.
# =============================================================================

library(tidyverse)
library(terra)
library(fuzzyjoin)
library(glue)
library(future)
library(furrr)
library(sf)
library(exactextractr)

plan(sequential)

source("scripts/helpers.R") # provides extract_years_from_files(), pull_focal_sf()

# Adjust these paths to your local or HPC setup before running.
home <- "H:/"
INPUT_DIR <- "H:/forestry-leakage/input"
OUTPUT_DIR <- "H:/forestry-leakage/output/run_202508"

# ── Enumerate treated-belt rasters and annual woody-vegetation layers ─────────
treated <- list.files(
  file.path(OUTPUT_DIR, "treated_rast"),
  full.names = TRUE
) %>%
  extract_years_from_files(
    pattern = "_([[:alnum:]]+)_\\d+\\_Y(\\d{4})",
    extension = ".tif"
  )

woody <- list.files(
  file.path(INPUT_DIR, "woody_veg/annual"),
  full.names = TRUE
) %>%
  extract_years_from_files(
    pattern = "([[:alnum:]]+)_(\\d{4})",
    extension = ".tif"
  ) %>%
  dplyr::select(filename, year) %>%
  rename(woody_filename = filename)

# NOTE (triage 1): source lists unmatch_data_buffer twice; neighbour omitted.
unmatch_data <- c(
  list.files(glue("{OUTPUT_DIR}/unmatch_data_buffer"), full.names = TRUE),
  list.files(glue("{OUTPUT_DIR}/unmatch_data_buffer"), full.names = TRUE)
) %>%
  extract_years_from_files(pattern = "_([[:alnum:]]+)_\\d+\\_Y(\\d{4})")

slats_name <- list.files(
  glue("{home}/nsw-habitat-loss/output/slats"),
  pattern = "*.tif$"
) %>%
  sapply(basename) %>%
  sapply(tools::file_path_sans_ext) %>%
  unname()
slats_rast_resampled <- rast(glue("{INPUT_DIR}/resampled/slats_resampled.tif"))
names(slats_rast_resampled) <- slats_name

plan(multisession, workers = 6)

# ── Match each treated raster to the nearest annual woody layer <= its year ──
joined_df <- treated %>%
  fuzzy_left_join(woody, by = c("year" = "year"), match_fun = list(`<=`)) %>%
  group_by(filename, year.x) %>%
  mutate(year.y = as.numeric(year.y)) %>%
  slice_max(order_by = -year.y, n = 1, with_ties = FALSE) %>%
  ungroup()
joined_df$run_id <- seq_len(nrow(joined_df))

unique_woody <- joined_df %>% distinct(year.y, .keep_all = TRUE)

# Reclassify the woody-veg raster to a 0/1 mask (classes 1,2,3 -> 1).
woody_calc <- future_map(seq_len(nrow(unique_woody)), function(i) {
  filename <- unique_woody$woody_filename[i]
  year <- unique_woody$year.y[i]
  woody_rast <- rast(filename)
  rcl <- matrix(c(0, 0, 1, 1, 2, 1, 3, 1), byrow = TRUE, ncol = 2)
  out_file <- as.character(glue(
    "{OUTPUT_DIR}/intermediate/woody_veg_class_{year}.tif"
  ))
  if (!file.exists(out_file)) {
    classify(woody_rast, rcl, others = 0, filename = out_file, overwrite = TRUE)
  }
  out_file
})

unique_woody_df <- unique_woody %>%
  dplyr::select(year.y) %>%
  mutate(woody_veg_class = woody_calc)
joined_df <- joined_df %>% left_join(unique_woody_df, by = "year.y")

# ── Focal state forest polygons ─────────────────────────────────────────────
focal_sf <- pull_focal_sf(unmatch_data)

focal_sf_df <- focal_sf %>%
  separate(run_id, c("run", "year")) %>%
  fuzzy_left_join(
    dplyr::select(joined_df, run, year.x, woody_veg_class),
    by = c("run", "year" = "year.x"),
    match_fun = list(`<=`)
  ) %>%
  rename(run = run.x) %>%
  group_by(run, year) %>%
  mutate(year.x = as.numeric(year.x)) %>%
  slice_max(order_by = -year.x, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(run, year, woody_veg_class)

focal_sf_tidy <- focal_sf %>%
  separate(run_id, c("run", "year")) %>%
  left_join(focal_sf_df)

# ── Woody vegetation inside each focal SF at the revocation year ─────────────
plan(multisession, workers = 8)
focal_sf_woody <- future_map(seq_len(nrow(focal_sf_tidy)), function(i) {
  row <- focal_sf_tidy[i, ]
  if ("geom" %in% names(row)) {
    st_geometry(row) <- "geom"
  }
  woody_rast <- rast(joined_df$woody_veg_class[[i]]) # NOTE (triage 2): positional
  area <- exactextractr::exact_extract(
    woody_rast,
    row,
    fun = "sum",
    max_cells_in_memory = 3e5
  ) *
    res(woody_rast)[1] *
    res(woody_rast)[2]
  sum(area)
}, .progress = TRUE)
plan(sequential)
focal_sf$woody_area_start_m2 <- unlist(focal_sf_woody)

# ── SLATS class-4 (woody) clearing inside each focal SF, per year ────────────
plan(multisession, workers = 4)
slats_resampled_wrap <- wrap(slats_rast_resampled)
focal_sf_slats <- future_map(seq_len(nrow(focal_sf_tidy)), function(i) {
  slats_rast_resampled <- unwrap(slats_resampled_wrap)
  row <- focal_sf_tidy[i, ] %>% st_cast("POLYGON")
  if ("geom" %in% names(row)) {
    st_geometry(row) <- "geom"
  }
  exact_extract(
    slats_rast_resampled,
    row,
    coverage_area = TRUE,
    fun = function(x, coverage_area) {
      bind_cols(x, coverage_area = coverage_area) %>%
        mutate_at(
          names(slats_rast_resampled),
          ~ if_else(is.na(.x) | (.x != 4), 0, (.x == 4) * coverage_area)
        ) %>%
        colSums(na.rm = TRUE)
    },
    force_df = TRUE,
    max_cells_in_memory = 3e6
  )
}, .progress = TRUE)
plan(sequential)

focal_sf_slats_sum <- focal_sf_slats %>%
  map(\(x) rowSums(x)) %>%
  reduce(cbind) %>%
  t() %>%
  as.data.frame()

deforestation_area <- focal_sf_slats_sum[, seq_along(names(slats_rast_resampled))]
names(deforestation_area) <- names(slats_rast_resampled)

# 1988-2005 SLATS layers are biannual; split evenly into annual values.
biannual_deforestation <- deforestation_area[, rep(1:9, each = 2)] / 2
names(biannual_deforestation) <- paste0("slats_", 1988:2005)
annual_deforestation <- deforestation_area[, 10:ncol(deforestation_area)]
names(annual_deforestation) <- paste0("slats_", 2006:2020)
ann_deforestation_area <- bind_cols(biannual_deforestation, annual_deforestation)
rownames(ann_deforestation_area) <- NULL

# On-site avoided clearing = clearing in the 5 years before revocation minus
# the 5 years after (revocation year y excluded from both).
five_year_deforestation <- c()
years <- as.numeric(joined_df$year.x) # NOTE (triage 2): positional alignment
for (i in seq_along(years)) {
  y <- years[i]
  pre <- sum(ann_deforestation_area[i, paste0("slats_", (y - 5):(y - 1))])
  post <- sum(ann_deforestation_area[i, paste0("slats_", (y + 1):(y + 5))])
  five_year_deforestation <- c(five_year_deforestation, pre - post)
}
focal_sf$five_year_clearing_ha <- unlist(five_year_deforestation * 0.0001)

st_write(
  focal_sf,
  file.path(OUTPUT_DIR, "focal_state_forest_gpkg/focal_sf_year.gpkg"),
  append = FALSE,
  delete_dsn = TRUE,
  delete_layer = TRUE
)

# ── Woody vegetation inside the treated leakage belt -> year_run_df.csv ──────
plan(multisession)
nsw_lga <- st_read(glue(
  "{INPUT_DIR}/nsw_border/NSW_LGA_GDA2020/NSW_LGA_GDA2020.shp"
))
nsw_border <- st_union(nsw_lga)

treated_polygon <- map(
  seq_len(nrow(joined_df)),
  safely(function(i) {
    treated_rast <- rast(joined_df$filename[[i]])
    out_name <- tools::file_path_sans_ext(basename(joined_df$filename[[i]]))
    out_file <- glue("{OUTPUT_DIR}/intermediate/{out_name}.shp")
    if (!file.exists(out_file)) {
      poly <- st_as_sf(as.polygons(treated_rast)) %>%
        st_transform(st_crs(nsw_border)) %>%
        st_intersection(nsw_border)
      st_write(poly, out_file)
    }
  }),
  .progress = TRUE
)
plan(sequential)

area <- map(
  seq_len(nrow(joined_df)),
  safely(function(i) {
    woody_rast <- rast(joined_df$woody_veg_class[[i]])
    out_name <- tools::file_path_sans_ext(basename(joined_df$filename[[i]]))
    out_file <- glue("{OUTPUT_DIR}/intermediate/{out_name}.shp")
    treated_polygon <- st_read(out_file) %>%
      filter(layer == 1) %>%
      summarise()
    exactextractr::exact_extract(
      woody_rast,
      treated_polygon,
      fun = "sum",
      max_cells_in_memory = 3e5
    ) *
      res(woody_rast)[1] *
      res(woody_rast)[2]
  }),
  .progress = TRUE
)

area_result <- lapply(area, \(x) sum(x$result))
write_rds(area, file.path(OUTPUT_DIR, "treated_rast/area.rds"))
joined_df$woody_veg_area <- unlist(area_result)

joined_df %>%
  dplyr::select(-year.y) %>%
  rename(year = year.x) %>%
  dplyr::select(-woody_veg_class) %>%
  write_csv(file.path(OUTPUT_DIR, "plot_data/year_run_df.csv"))

plan(sequential)
