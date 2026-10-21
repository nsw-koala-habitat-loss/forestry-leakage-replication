# =============================================================================
# Script 08a: Export the inputs for the matching-design map (match_plots_small)
# =============================================================================
# [Reference] Run ONCE by the authors on the machine that holds the full
# analysis outputs (H:/forestry-leakage). Readers do not need to run this: its
# outputs ship in data/match_map/ and are read by Script 08.
#
# It writes small, cropped copies of every layer the map draws, taken from the
# same objects that scripts/plots.R used to make plots/match_plots_small.png:
#
#   data/match_map/match_map_layers.gpkg
#     plots_buffer       Matched grid cells, Pilliga East (panel A)
#     plots_neighbour    Matched grid cells, Coolangubra (panel B)
#     state_forests      State Forest polygons (site_id, SFName, revocation_year)
#     protected_areas    NPWS estate, filtered by gazettal year as in plots.R
#   data/match_map/woody_veg_buffer.tif     Woody vegetation, panel A extent
#   data/match_map/woody_veg_neighbour.tif  Woody vegetation, panel B extent
#   data/match_map/match_map_panels.csv     Site IDs and names for each panel
# =============================================================================

library(glue)
library(sf)
library(dplyr)
library(terra)
library(readr)

INPUT_DIR <- file.path("H:/forestry-leakage/input")
OUTPUT_DIR <- file.path("H:/forestry-leakage/output/run_202508")
REP_DIR <- file.path("H:/forestry-leakage/forestry_leakage_rep") # set to the package root

OUT_DIR <- file.path(REP_DIR, "data", "match_map")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
gpkg_path <- file.path(OUT_DIR, "match_map_layers.gpkg")
if (file.exists(gpkg_path)) file.remove(gpkg_path)

# Same objects as scripts/plots.R
if (!file.exists(glue("{OUTPUT_DIR}/plot_data/plot_data.RData"))) {
  source(glue("{INPUT_DIR}/../scripts/load_plot_data.R"))
} else {
  load(glue("{OUTPUT_DIR}/plot_data/plot_data.RData"))
}
wv <- rast(glue("{INPUT_DIR}/woody_veg/woody_nsw.tif"))

panels <- tibble(
  panel = c("A", "B"),
  channel = c("buffer", "neighbour"),
  SFName = c(rep_buffer, rep_neighbour),
  site_id = c(rep_buffer_id, rep_neighbour_id)
)
match_data_sp <- list(
  buffer = match_data_sp_buffer,
  neighbour = match_data_sp_neighbour
)

# Map extent: focal forest + 20 km, squared (as in plots.R)
panel_bbox <- function(id) {
  state_forests_revoked %>%
    filter(site_id == id) %>%
    st_union() %>%
    st_buffer(20000) %>%
    get_square_bbox()
}
bboxes <- lapply(panels$site_id, panel_bbox)
names(bboxes) <- panels$channel
bbox_union <- bboxes %>%
  lapply(\(b) st_as_sfc(b)) %>%
  do.call(c, .) %>%
  st_union()

# Matched plots: full matched sample for each representative site
for (ch in panels$channel) {
  st_write(match_data_sp[[ch]], gpkg_path, layer = glue("plots_{ch}"), quiet = TRUE)
}

# State Forests: whole polygons touching either panel
state_forests_classed %>%
  select(site_id, SFName, revocation_year) %>%
  st_intersects_subset(bbox_union) %>%
  st_write(gpkg_path, layer = "state_forests", quiet = TRUE)

# Protected areas: treatment_control_plot() filters with
#   filter(GAZ_DATE < as.numeric(revocation_year) + 5)
# where `revocation_year` is not a column of protected_areas, so it resolves to
# whatever `revocation_year` exists in the session. Apply the identical
# expression here and record the value so the map is reproduced exactly.
stopifnot(!"revocation_year" %in% names(protected_areas))
if (!exists("revocation_year", envir = globalenv())) {
  stop("`revocation_year` not found in the session: the protected-area filter in ",
       "treatment_control_plot() cannot be reproduced. Check how plots.R was run.")
}
pa_cutoff_year <- as.numeric(get("revocation_year", envir = globalenv()))
message("Protected-area filter value(s) of revocation_year: ",
        paste(unique(pa_cutoff_year), collapse = ", "))

protected_areas %>%
  mutate(GAZ_DATE = substr(GAZ_DATE, 0, 4) %>% as.numeric()) %>%
  filter(GAZ_DATE < pa_cutoff_year + 5) %>%
  st_intersects_subset(bbox_union) %>%
  st_write(gpkg_path, layer = "protected_areas", quiet = TRUE)

# Woody vegetation: crop to each panel extent exactly as treatment_control_plot()
for (ch in panels$channel) {
  bbox_terra <- bboxes[[ch]] %>%
    st_as_sfc() %>%
    vect() %>%
    project(wv) %>%
    ext()
  writeRaster(
    crop(wv, bbox_terra),
    file.path(OUT_DIR, glue("woody_veg_{ch}.tif")),
    overwrite = TRUE,
    gdal = c("COMPRESS=DEFLATE")
  )
}

panels %>%
  mutate(pa_cutoff_year = paste(unique(pa_cutoff_year), collapse = ";")) %>%
  write_csv(file.path(OUT_DIR, "match_map_panels.csv"))

# Summary for the log
print(st_layers(gpkg_path))
print(list.files(OUT_DIR, full.names = TRUE) %>% file.size() %>% setNames(list.files(OUT_DIR)))
