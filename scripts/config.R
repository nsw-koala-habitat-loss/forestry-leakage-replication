# Replication package configuration.
# All paths are relative to the forestry_leakage_rep/ root directory.
# Run scripts from that directory (e.g. setwd("forestry_leakage_rep")).

BASE_DIR  <- getwd()
DATA_DIR  <- file.path(BASE_DIR, "data")
PLOT_DIR  <- file.path(BASE_DIR, "plots")

if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR)
