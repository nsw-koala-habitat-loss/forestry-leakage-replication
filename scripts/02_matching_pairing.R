# =============================================================================
# Script 02: Treatment–control pairing via propensity-score matching
# =============================================================================
# PURPOSE
#   For each revoked state forest (focal SF), constructs a treatment group
#   (grid cells within the leakage zone) and a donor control pool, then runs
#   19 MatchIt specifications. Outputs paired match data used in Script 03.
#
#   Two leakage channels are modelled:
#     buffer   – leakage to private native forestry within 10 km of the focal SF
#     neighbour – leakage to neighbouring state forests within 10 km
#
# REFERENCE SCRIPT
#   Requires the full spatial inputs and the covariate Parquet files produced
#   by Script 01. Run on HPC; not directly executable with this package.
#
# INPUTS (HPC paths)
#   {INPUT_DIR}/state_forests/state_forests_revoked_update.gpkg
#   {INPUT_DIR}/nsw_grid/nsw_grid_points.parquet
#   {INPUT_DIR}/matching_data/parquet/*.parquet        (from Script 01)
#   {INPUT_DIR}/protected_areas/…/NPWS_EstateInternalBoundaries.shp
#
# OUTPUTS  (per focal SF, in OUTPUT_DIR)
#   match_data_pair_{buffer,neighbour}/match_data_long_all_{i}_{name}.rds
#   matchit_{buffer,neighbour}/matchit_{buffer,neighbour}_summary_{i}_{name}.csv
#   match_data_pair_{buffer,neighbour}/match_data_long_{i}_{name}.parquet
# =============================================================================

library(tidyverse)
library(sf)
library(terra)
library(arrow)
library(glue)
library(MatchIt)
library(INLA)
library(purrr)
library(geoarrow)
library(car)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && !is.na(args[1])) {
  job_id <- as.numeric(args[1])
}

if (exists("job_id")) {
  home <- "/mnt/hpccs01/home/<your-hpc-username>"
} else {
  home <- "H:" # adjust to your local path
}

INPUT_DIR <- file.path(home, "forestry-leakage/input")
OUTPUT_DIR <- file.path(home, "forestry-leakage/output/run_202508")

# -- Matching formulae --------------------------------------------------------
buffer_formula <- treat ~ DistCity +
  DistRoad +
  EcolCond +
  Remote_3 +
  LandUse_Natural +
  LandUse_Dryland +
  LandUse_Irrigated +
  LandUse_Intensive +
  NatVegReg_Cat1 +
  NatVegReg_Cat2 +
  Eucalypt2018 +
  Rainforest2018 +
  Plantation2018 +
  temp_q2 +
  temp_q3 +
  temp_q4 +
  prec_q2 +
  prec_q3 +
  prec_q4 +
  slope +
  ScEc_PC1 +
  ScEc_PC2 +
  Soil_PC1 +
  Soil_PC2 +
  Fire_1 +
  Fire_2 +
  Drought

neighbour_formula <- treat ~ DistCity +
  DistRoad +
  EcolCond +
  Remote_3 +
  Eucalypt2018 +
  Rainforest2018 +
  Plantation2018 +
  temp_q2 +
  temp_q3 +
  temp_q4 +
  prec_q2 +
  prec_q3 +
  prec_q4 +
  slope +
  ScEc_PC1 +
  ScEc_PC2 +
  Soil_PC1 +
  Soil_PC2 +
  Fire_1 +
  Fire_2 +
  Drought

# -- 19 matching specifications to evaluate -----------------------------------
spec_list <- list(
  list(method = "nearest", distance = "glm", link = "logit"),
  list(method = "nearest", distance = "glm", link = "probit"),
  list(method = "nearest", distance = "robust_mahalanobis"),
  list(method = "nearest", caliper = 0.9),
  list(method = "nearest", caliper = 0.5),
  list(method = "nearest", caliper = 0.1),
  list(method = "nearest", replace = TRUE),

  list(
    method = "nearest",
    distance = "glm",
    link = "logit",
    new_vars = c("WoodyVegTreat")
  ),
  list(
    method = "nearest",
    distance = "glm",
    link = "probit",
    new_vars = c("WoodyVegTreat")
  ),
  list(
    method = "nearest",
    distance = "robust_mahalanobis",
    new_vars = c("WoodyVegTreat")
  ),
  list(method = "nearest", caliper = 0.9, new_vars = c("WoodyVegTreat")),
  list(method = "nearest", caliper = 0.5, new_vars = c("WoodyVegTreat")),
  list(method = "nearest", caliper = 0.1, new_vars = c("WoodyVegTreat")),

  list(
    method = "nearest",
    distance = "glm",
    link = "logit",
    new_vars = c("WoodyVegTreat", "PreClearing")
  ),
  list(
    method = "nearest",
    distance = "glm",
    link = "probit",
    new_vars = c("WoodyVegTreat", "PreClearing")
  ),
  list(
    method = "nearest",
    distance = "robust_mahalanobis",
    new_vars = c("WoodyVegTreat", "PreClearing")
  ),
  list(
    method = "nearest",
    caliper = 0.9,
    new_vars = c("WoodyVegTreat", "PreClearing")
  ),
  list(
    method = "nearest",
    caliper = 0.5,
    new_vars = c("WoodyVegTreat", "PreClearing")
  ),
  list(
    method = "nearest",
    caliper = 0.1,
    new_vars = c("WoodyVegTreat", "PreClearing")
  )
)

# -- Helper: create output subdirectories -------------------------------------
create_output_dir <- function(parent_dir, sub_dir_list) {
  for (sub_dir in sub_dir_list) {
    full_path <- file.path(parent_dir, sub_dir)
    if (!dir.exists(full_path)) dir.create(full_path, recursive = TRUE)
  }
}

# -- Helper: reclassify covariates into binary/quartile indicators ------------
# Cutpoints are NSW-wide percentiles computed from the full 10 km grid (not the
# matched sample), so all focal SFs are matched on a common scale.
add_categorical_vars <- function(df) {
  temp_25 <- 16.967110
  temp_50 <- 18.684704
  temp_75 <- 20.249214
  prec_25 <- 233.15067
  prec_50 <- 369.97745
  prec_75 <- 594.03262

  df %>%
    mutate(
      Remote_3 = if_else(remote2016 >= 3, 1, 0),
      LandUse_Natural = if_else(LandUse == 1, 1, 0),
      LandUse_Dryland = if_else(LandUse == 2, 1, 0),
      LandUse_Irrigated = if_else(LandUse == 3, 1, 0),
      LandUse_Intensive = if_else(LandUse == 4, 1, 0),
      NatVegReg_Cat1 = if_else(NatVegReg == 1, 1, 0),
      NatVegReg_Cat2 = if_else(NatVegReg %in% c(2, 3), 1, 0),
      Fire_1 = if_else(is.na(Fire), 0, if_else(Fire == 1, 1, 0)),
      Fire_2 = if_else(is.na(Fire), 0, if_else(Fire == 2, 1, 0))
    ) %>%
    mutate(
      temp_q1 = if_else(temp < temp_25, 1, 0),
      temp_q2 = if_else(between(temp, temp_25, temp_50), 1, 0),
      temp_q3 = if_else(between(temp, temp_50, temp_75), 1, 0),
      temp_q4 = if_else(temp > temp_75, 1, 0),
      prec_q1 = if_else(prec < prec_25, 1, 0),
      prec_q2 = if_else(between(prec, prec_25, prec_50), 1, 0),
      prec_q3 = if_else(between(prec, prec_50, prec_75), 1, 0),
      prec_q4 = if_else(prec > prec_75, 1, 0)
    )
}

# -- Helper: stepwise VIF removal with collinearity pre-check -----------------
find_perfectly_collinear_pairs <- function(cor_matrix, cor_threshold) {
  upper_tri <- upper.tri(cor_matrix, diag = FALSE)
  high_cor_pairs <- which(
    abs(cor_matrix) >= cor_threshold & upper_tri,
    arr.ind = TRUE
  )
  vars_to_remove <- character()
  if (nrow(high_cor_pairs) > 0) {
    for (i in seq_len(nrow(high_cor_pairs))) {
      vars_to_remove <- c(
        vars_to_remove,
        colnames(cor_matrix)[high_cor_pairs[i, 2]]
      )
    }
  }
  unique(vars_to_remove)
}

stepwise_vif_with_cor_check <- function(
  formula,
  data,
  vif_threshold = 5,
  cor_threshold = 0.6
) {
  formula_vars <- attr(terms(formula), "term.labels")
  response_var <- as.character(formula)[2]
  constant_vars <- formula_vars[sapply(formula_vars, function(v) {
    length(unique(data[[v]])) <= 1
  })]

  if (length(constant_vars) > 0) {
    message(
      "Removing constant variables: ",
      paste(constant_vars, collapse = ", ")
    )
    formula_vars <- setdiff(formula_vars, constant_vars)
  }

  predictors <- formula_vars

  if (length(predictors) >= 2) {
    numeric_predictors <- predictors[sapply(data[, predictors], is.numeric)]
    if (length(numeric_predictors) >= 2) {
      cor_matrix <- cor(data[, numeric_predictors])
      diag(cor_matrix) <- 0
      collinear <- find_perfectly_collinear_pairs(cor_matrix, cor_threshold)
      if (length(collinear) > 0) {
        message(
          "Removing collinear variables: ",
          paste(collinear, collapse = ", ")
        )
        predictors <- setdiff(predictors, collinear)
      }
    }
  }

  while (TRUE) {
    if (length(predictors) < 2) {
      final_formula <- as.formula(paste(
        response_var,
        "~",
        paste(predictors, collapse = " + ")
      ))
      final_model <- glm(final_formula, data = data, family = "binomial")
      return(list(model = final_model, formula = final_formula))
    }
    current_formula <- as.formula(paste(
      response_var,
      "~",
      paste(predictors, collapse = " + ")
    ))
    model <- glm(current_formula, data = data, family = "binomial")
    vifs <- vif(model)
    max_vif <- max(vifs)
    if (max_vif < vif_threshold) {
      message("All VIFs below threshold ", vif_threshold)
      return(list(model = model, formula = current_formula))
    }
    message(
      "Removing '",
      names(which.max(vifs)),
      "' (VIF = ",
      round(max_vif, 2),
      ")"
    )
    predictors <- predictors[predictors != names(which.max(vifs))]
  }
}

# -- Helper: evaluate whether a MatchIt object meets balance criteria ---------
matchit_meets_criteria <- function(
  obj,
  parallel_trends_p = 0.05,
  SMD_threshold = 0.25,
  prop_to_keep = 0.8
) {
  if (is.null(obj)) {
    return(c(FALSE, FALSE, FALSE))
  }

  parallel_trends_test <- function(obj) {
    x <- obj$match_data_long %>%
      filter(t < 0) %>%
      group_by(plot_id) %>%
      arrange(plot_id, t) %>%
      mutate(cum_prop_clearing = cumsum(prop_clearing)) %>%
      ungroup() %>%
      inla(
        cum_prop_clearing ~ 1 + t + CI + CI * t + f(plot_id, model = "iid"),
        data = .,
        family = "gaussian",
        num.threads = 1
      )
    # Two-tailed Bayesian equivalence check: the test PASSES (returns TRUE) if
    # the credible interval of the treatment-trend interaction (t:CI) straddles
    # zero — i.e., we cannot rule out parallel trends. It FAILS if the CI lies
    # entirely above or entirely below zero, indicating diverging pre-trends.
    sig_pos <- (inla.qmarginal(parallel_trends_p, x$marginals.fixed$`t:CI`) >
      0) &
      (inla.qmarginal(0.5, x$marginals.fixed$`t:CI`) > 0)
    sig_neg <- (inla.qmarginal(0.5, x$marginals.fixed$`t:CI`) < 0) &
      (inla.qmarginal(1 - parallel_trends_p, x$marginals.fixed$`t:CI`) < 0)
    !(sig_pos | sig_neg)
  }

  parallel_trends_pass <- purrr::possibly(
    parallel_trends_test,
    otherwise = FALSE
  )(obj)

  matchit_summary <- summary(obj$matchit_obj, std.dev = "pooled")
  smd <- matchit_summary$sum.matched[, "Std. Mean Diff."]
  smd[is.na(smd)] <- 0
  smd_pass <- all(abs(smd) < SMD_threshold)
  prop_pass <- (matchit_summary$nn["Matched", "Treated"] /
    matchit_summary$nn["All", "Treated"]) >
    prop_to_keep

  c(parallel_trends_pass, smd_pass, prop_pass)
}

# -- Helper: run MatchIt and return tidy match data ---------------------------
matchit_param <- function(unmatch_data_obj, formula, ...) {
  args <- as.list(match.call())
  if ("new_vars" %in% names(args)) {
    formula <- update(formula, reformulate(c(".", args$new_vars)))
  }

  unmatch_data_sp <- unmatch_data_obj$unmatch_data
  focal_rev_year <- unmatch_data_obj$focal_rev_year
  woody_veg_lyr_i <- unmatch_data_obj$woody_veg_lyr

  unmatch_data <- st_drop_geometry(unmatch_data_sp) %>%
    drop_na(all_of(all.vars(formula)))

  m_out1 <- matchit(formula, data = unmatch_data, ...)
  match_data <- match_data(m_out1, subclass = "subclass")

  if (!("subclass" %in% colnames(match_data))) {
    match_data <- mutate(match_data, subclass = NA)
  }

  match_data_sp <- unmatch_data_sp %>%
    inner_join(match_data, by = "plot_id", suffix = c("", ".y")) %>%
    dplyr::select(-dplyr::ends_with(".y"))

  match_data_long <- match_data %>%
    pivot_longer(
      starts_with("slats_"),
      names_prefix = "slats_",
      names_to = "year",
      values_to = "forest_clearing"
    ) %>%
    mutate(year = as.numeric(year)) %>%
    filter(between(year, focal_rev_year - 5, focal_rev_year + 5)) %>%
    mutate(
      t = year - focal_rev_year,
      # focal_rev_year itself is dropped just below (filter(t != 0)); `>` keeps
      # the before/after split unambiguous if that filter is ever relaxed.
      BA = if_else(year > focal_rev_year, 1, 0),
      CI = if_else(treat == 1, 1, 0)
    ) %>%
    filter(t != 0) %>%
    mutate(
      prop_clearing = if_else(
        !!sym(woody_veg_lyr_i) == 0 | !!sym(woody_veg_lyr_i) < forest_clearing,
        0,
        forest_clearing / !!sym(woody_veg_lyr_i)
      )
    )

  list(
    match_data = match_data,
    match_data_long = match_data_long,
    matchit_obj = m_out1,
    match_data_sp = match_data_sp
  )
}

# -- Initialise shared objects (loads lazily on first call) -------------------
# <<- assignments are intentional: init() populates the calling environment so
# that helpers like st_intersects_subset can reference shared spatial objects
# without requiring them to be passed as arguments to every match_* call.
small_init <- function() {
  INPUT_DIR <<- file.path(home, "forestry-leakage/input")
  OUTPUT_DIR <<- file.path(home, "forestry-leakage/output/run_202508")
  if (!exists("state_forests_revoked")) {
    state_forests_revoked <<- st_read(
      file.path(INPUT_DIR, "state_forests/state_forests_revoked_update.gpkg"),
      quiet = TRUE
    )
  }
  state_forests_revoked_subset <<- filter(state_forests_revoked, diff < 0)
}

init <- function(run_spatial = TRUE) {
  buffer_area_size <<- 10000

  INPUT_DIR <<- file.path(home, "forestry-leakage/input")
  OUTPUT_DIR <<- file.path(home, "forestry-leakage/output/run_202508")

  create_output_dir(
    OUTPUT_DIR,
    c(
      "focal_state_forest_gpkg",
      "match_data_pair_buffer",
      "match_data_pair_neighbour",
      "match_data_sp_buffer",
      "match_data_sp_neighbour",
      "matchit_buffer",
      "matchit_neighbour",
      "plot_data",
      "plots",
      "treated_polygon_gpkg",
      "unmatch_data_buffer",
      "unmatch_data_neighbour"
    )
  )

  matching_df <<- open_dataset(file.path(INPUT_DIR, "matching_data/parquet"))

  if (!exists("state_forests_revoked")) {
    state_forests_revoked <<- st_read(
      file.path(INPUT_DIR, "state_forests/state_forests_revoked_update.gpkg"),
      quiet = TRUE
    )
  }
  state_forests_revoked_subset <<- filter(state_forests_revoked, diff < 0)

  if (run_spatial) {
    if (!exists("state_forests_all")) {
      state_forests_all <<- read_sf(
        file.path(INPUT_DIR, "state_forests/state_forests_all.shp"),
        quiet = TRUE
      )
    }

    if (!exists("state_forests_clean")) {
      state_forests_clean <<- state_forests_all %>%
        dplyr::rename(
          CurrentTen = CrrntTn,
          Revocation = Revoctn,
          revocation_year = rvctn_y,
          area_km2 = are_km2
        ) %>%
        dplyr::select(
          SFName,
          SFNo,
          CurrentTen,
          Revocation,
          revoked,
          revocation_year,
          area_km2
        )
    }

    if (!exists("grid_points")) {
      grid_points <<- open_dataset(file.path(
        INPUT_DIR,
        "nsw_grid/nsw_grid_points.parquet"
      ))
    }

    if (!exists("protected_areas")) {
      protected_areas <<- st_read(
        file.path(
          INPUT_DIR,
          "protected_areas/npws_estateinternalboundaries/NPWS_EstateInternalBoundaries.shp"
        ),
        quiet = TRUE
      ) %>%
        st_transform(st_crs(state_forests_revoked)) %>%
        mutate(GAZ_YEAR = as.numeric(substr(GAZ_DATE, 1, 4)))
    }
  }

  st_intersects_subset <<- function(a, b, negate = FALSE, within = FALSE) {
    intersect_list <- if (within) st_within(a, b) else st_intersects(a, b)
    intersects <- lengths(intersect_list) > 0
    if (negate) a[!intersects, ] else a[intersects, ]
  }

  if (run_spatial) {
    if (!exists("state_forest_mask")) {
      state_forest_mask <<- summarise(state_forests_clean)
    }
  }

  matchit_balance_criteria <<- function(matchit_list) {
    purrr::map(
      matchit_list,
      matchit_meets_criteria,
      .progress = "MatchIt balance evaluation"
    ) %>%
      do.call(rbind, .)
  }

  matchit_list_summary <<- function(matchit_list, spec_list) {
    spec_df <- spec_list %>%
      lapply(function(x) {
        x$new_vars <- glue_collapse(x$new_vars, sep = ", ")
        as.data.frame(x)
      }) %>%
      bind_rows() %>%
      mutate(
        spec_id = as.character(row_number()),
        distance = if_else(is.na(distance), "glm", distance),
        link = if_else(is.na(link), "logit", link),
        replace = if_else(is.na(replace), FALSE, replace)
      )

    balance_criteria <- matchit_balance_criteria(matchit_list)
    colnames(balance_criteria) <- c("ParallelTrends", "SMD", "PropKeep")
    matchit_idx <- which(
      rowSums(balance_criteria) == max(rowSums(balance_criteria))
    )[1]

    matchit_summary <- lapply(matchit_list, \(x) {
      if (is.null(x)) {
        return(NA)
      }
      s <- summary(x$matchit_obj)$sum.matched %>% as.data.frame()
      s$covariate <- row.names(s)
      row.names(s) <- NULL
      s
    })
    names(matchit_summary) <- seq_along(matchit_summary)

    balance_criteria_df <- as.data.frame(balance_criteria) %>%
      mutate(spec_id = as.character(row_number()))

    bind_rows(matchit_summary, .id = "spec_id") %>%
      left_join(spec_df, by = "spec_id") %>%
      mutate(model_selected = as.numeric(spec_id) == matchit_idx) %>%
      left_join(balance_criteria_df, by = "spec_id")
  }
}

# -- Main function: buffer leakage matching -----------------------------------
# Treatment group: grid cells on private land within 10 km of the focal SF that
# were not yet protected at the time of revocation.
# Donor pool: cells outside all state forests, not near any recently-revoked SF,
# and within 500 km — private-land comparators only.
match_buffer <- function(
  i,
  run_spatial = FALSE,
  write_gpkg = TRUE,
  evaluate_matchit = FALSE
) {
  small_init()

  n_individual_sf <- nrow(state_forests_revoked_subset)
  if (i <= n_individual_sf) {
    focal_sf <- state_forests_revoked_subset[i, ]
    focal_rev_year <- focal_sf[["revocation_year"]]
    focal_sf_name <- focal_sf[["SFName"]]
  } else {
    unique_years <- sort(unique(state_forests_revoked$revocation_year))
    focal_rev_year <- unique_years[i - n_individual_sf]
    focal_sf <- state_forests_revoked %>%
      filter(as.numeric(revocation_year) == focal_rev_year) %>%
      summarise()
    focal_sf_name <- glue("Y{focal_rev_year}")
  }

  print(glue("Buffer Iteration {i}: {focal_sf_name}"))

  formula <- buffer_formula
  unmatch_data_path <- glue(
    "{OUTPUT_DIR}/unmatch_data_buffer/unmatch_data_buffer_{i}_{focal_sf_name}.rds"
  )

  if (run_spatial | !file.exists(unmatch_data_path)) {
    init(TRUE)

    start_year <- focal_rev_year - 5
    woody_veg_names <- names(dplyr::select(
      matching_df,
      starts_with("woody_veg_")
    ))
    woody_veg_years <- as.numeric(str_split_i(woody_veg_names, "_", 3))

    # Select the annual woody-veg layer closest to start_year from *below* by
    # setting future-year differences to 99. This avoids using post-treatment
    # data when computing the pre-treatment woody-veg offset for the model.
    year_diff <- woody_veg_years - start_year
    year_diff[year_diff > 0] <- 99
    woody_veg_lyr_i <- woody_veg_names[which.min(abs(year_diff))]
    woody_veg_lyr_j <- woody_veg_names[which.min(abs(
      woody_veg_years - focal_rev_year
    ))]

    PA_names <- names(dplyr::select(matching_df, starts_with("PA_")))
    PA_years <- as.numeric(str_split_i(PA_names, "_", 2))
    PA_lyr_i <- PA_names[match(start_year, PA_years)] %||% "PA_1988"
    # Screen on the PA layer covering the end of the BACI panel (focal_rev_year
    # + 5). The annual PA layers are cumulative (protected_areas.R: GAZ_YEAR <
    # year, never de-gazetted), so PA_{y} flags every cell gazetted on/before
    # y - 1; PA_{focal_rev_year + 6} drops any cell protected through the panel end.
    PA_lyr_j <- PA_names[match(focal_rev_year + 6, PA_years)] %||% "PA_2016"

    matching_df_woody_veg <- matching_df %>%
      filter(!!sym(woody_veg_lyr_i) > 0 & !!sym(PA_lyr_j) < 0.5) %>%
      dplyr::select(plot_id) %>%
      collect()

    grid_points_i <- grid_points %>%
      filter(plot_id %in% matching_df_woody_veg$plot_id) %>%
      st_as_sf() %>%
      inner_join(matching_df_woody_veg, by = "plot_id")

    focal_sf_buffer <- st_buffer(focal_sf, buffer_area_size) %>%
      st_difference(focal_sf)
    current_pa <- protected_areas %>%
      filter(GAZ_DATE <= (focal_rev_year + 5)) %>%
      summarise()
    focal_treated <- focal_sf_buffer %>%
      st_difference(state_forest_mask) %>%
      summarise() %>%
      st_difference(current_pa) %>%
      summarise() %>%
      st_buffer(0)

    focal_grid_points <- st_intersects_subset(
      grid_points_i,
      focal_treated,
      within = TRUE
    )

    # Donor pool spatial filters (each line removes a category of invalid comparators):
    #   1. Must be outside all state forests (private land only for buffer channel).
    #   2. Must be outside the 10 km buffer of any SF revoked within the past 5 years
    #      (prevents contamination from concurrent leakage events).
    #   3. Must be outside the focal SF's own 10 km buffer (not a treated unit).
    #   4. Must be within 500 km of the focal SF (broad geographic comparability).
    donor_pool_points <- grid_points_i %>%
      st_intersects_subset(state_forests_clean, negate = TRUE) %>%
      st_intersects_subset(
        state_forests_clean %>%
          filter(between(
            revocation_year,
            focal_rev_year - 5,
            focal_rev_year
          )) %>%
          st_buffer(buffer_area_size),
        negate = TRUE
      ) %>%
      st_intersects_subset(
        st_buffer(focal_sf, buffer_area_size),
        negate = TRUE
      ) %>%
      st_intersects_subset(st_buffer(focal_sf, 500000), within = TRUE)

    if (nrow(donor_pool_points) == 0) {
      return()
    }

    matching_df_subset <- matching_df %>%
      filter(
        plot_id %in% c(focal_grid_points$plot_id, donor_pool_points$plot_id)
      ) %>%
      mutate(
        WoodyVeg = !!sym(woody_veg_lyr_i),
        WoodyVegTreat = !!sym(woody_veg_lyr_j),
        WoodyVegDelta = !!sym(woody_veg_lyr_i) - !!sym(woody_veg_lyr_j),
        PA_Current = !!sym(PA_lyr_i),
        PA_Final = !!sym(PA_lyr_j)
      ) %>%
      collect()

    treatment <- focal_grid_points %>%
      left_join(matching_df_subset, by = "plot_id") %>%
      mutate(
        total_clearing = rowSums(across(starts_with("slats_"))),
        PreClearing = rowSums(across(any_of(paste0(
          "slats_",
          seq(focal_rev_year - 5, focal_rev_year - 1)
        ))))
      )

    control <- donor_pool_points %>%
      left_join(matching_df_subset, by = "plot_id") %>%
      mutate(
        total_clearing = rowSums(across(starts_with("slats_"))),
        PreClearing = rowSums(across(any_of(paste0(
          "slats_",
          seq(focal_rev_year - 5, focal_rev_year - 1)
        ))))
      )

    matching_data_sf <- bind_rows(
      treatment = treatment,
      control = control,
      .id = "treat"
    ) %>%
      mutate(
        treat = treat == "treatment",
        woody_veg_start = !!sym(woody_veg_lyr_i)
      ) %>%
      filter(PA_Final < 0.5) %>%
      drop_na(treat, woody_veg_start, any_of(all.vars(formula)))

    matching_data_out <- list(
      focal_sf = focal_sf,
      focal_sf_buffer = focal_sf_buffer,
      focal_rev_year = focal_rev_year,
      unmatch_data = matching_data_sf,
      woody_veg_lyr = woody_veg_lyr_i
    )
    write_rds(matching_data_out, unmatch_data_path, compress = "none")
  } else {
    init(FALSE)
    matching_data_out <- read_rds(unmatch_data_path)
  }

  matching_data_out$unmatch_data <- add_categorical_vars(
    matching_data_out$unmatch_data
  )

  no_collinear <- stepwise_vif_with_cor_check(
    formula,
    st_drop_geometry(matching_data_out$unmatch_data)
  )
  safe_match <- purrr::possibly(
    function(spec) {
      do.call(
        function(...) {
          matchit_param(matching_data_out, no_collinear$formula, ...)
        },
        spec
      )
    },
    otherwise = NULL
  )
  matchit_obj <- map(spec_list, safe_match, .progress = TRUE)

  match_data_long_path <- glue(
    "{OUTPUT_DIR}/match_data_pair_buffer/match_data_long_all_{i}_{focal_sf_name}.rds"
  )
  write_rds(
    lapply(matchit_obj, \(x) {
      dplyr::select(
        x$match_data_long,
        any_of(c(
          "plot_id",
          "treat",
          "year",
          "BA",
          "CI",
          "t",
          "forest_clearing",
          "prop_clearing",
          "subclass",
          "woody_veg_start",
          "KoalaHabitatSuitability"
        ))
      )
    }),
    match_data_long_path
  )

  if (!evaluate_matchit) {
    return()
  }

  matchit_summary <- matchit_list_summary(matchit_obj, spec_list)
  matchit_criteria_df <- matchit_summary %>%
    dplyr::select(any_of(c("spec_id", "ParallelTrends", "SMD", "PropKeep"))) %>%
    group_by(spec_id) %>%
    mutate(spec_id = as.numeric(spec_id)) %>%
    arrange(spec_id) %>%
    summarise(across(ParallelTrends:PropKeep, max)) %>%
    rowwise() %>%
    mutate(balance_criteria = sum(c(ParallelTrends, SMD, PropKeep)))

  matchit_criteria_idx <- matchit_criteria_df[
    matchit_criteria_df$balance_criteria ==
      max(matchit_criteria_df$balance_criteria),
    "spec_id"
  ][[1]][1]

  write_csv(
    matchit_summary,
    glue(
      "{OUTPUT_DIR}/matchit_buffer/matchit_buffer_summary_{i}_{focal_sf_name}.csv"
    )
  )
  write_parquet(
    matchit_obj[[matchit_criteria_idx]]$match_data_long,
    glue(
      "{OUTPUT_DIR}/match_data_pair_buffer/match_data_long_{i}_{focal_sf_name}.parquet"
    )
  )
  if (write_gpkg) {
    st_write(
      matchit_obj[[matchit_criteria_idx]]$match_data_sp,
      glue(
        "{OUTPUT_DIR}/match_data_sp_buffer/match_data_sp_{i}_{focal_sf_name}.gpkg"
      ),
      append = FALSE
    )
  }

  print(glue("Finished iteration {i}: {focal_sf_name}"))
}

# -- Main function: neighbour leakage matching --------------------------------
# Treatment group: grid cells inside neighbouring (intact) state forests within
# 10 km of the focal SF — leakage via increased logging pressure on adjacent SFs.
# Donor pool: cells inside intact state forests that are not near any recently-
# revoked SF and not within the focal SF's own buffer.
match_neighbour <- function(
  i,
  run_spatial = FALSE,
  write_gpkg = TRUE,
  evaluate_matchit = FALSE
) {
  small_init()

  n_individual_sf <- nrow(state_forests_revoked_subset)
  if (i <= n_individual_sf) {
    focal_sf <- state_forests_revoked_subset[i, ]
    focal_rev_year <- focal_sf[["revocation_year"]]
    focal_sf_name <- focal_sf[["SFName"]]
  } else {
    unique_years <- sort(unique(state_forests_revoked$revocation_year))
    focal_rev_year <- unique_years[i - n_individual_sf]
    focal_sf <- state_forests_revoked %>%
      filter(as.numeric(revocation_year) == focal_rev_year) %>%
      summarise()
    focal_sf_name <- glue("Y{focal_rev_year}")
  }

  print(glue("Neighbour Iteration {i}: {focal_sf_name}"))

  formula <- neighbour_formula
  unmatch_data_path <- glue(
    "{OUTPUT_DIR}/unmatch_data_neighbour/unmatch_data_neighbour_{i}_{focal_sf_name}.rds"
  )

  if (run_spatial | !file.exists(unmatch_data_path)) {
    init(TRUE)

    start_year <- focal_rev_year - 5
    woody_veg_names <- names(dplyr::select(
      matching_df,
      starts_with("woody_veg_")
    ))
    woody_veg_years <- as.numeric(str_split_i(woody_veg_names, "_", 3))

    # As in match_buffer: select the closest woody-veg layer from below start_year
    # to avoid using post-treatment data in the pre-treatment offset.
    year_diff <- woody_veg_years - start_year
    year_diff[year_diff > 0] <- 99
    woody_veg_lyr_i <- woody_veg_names[which.min(abs(year_diff))]
    woody_veg_lyr_j <- woody_veg_names[which.min(abs(
      woody_veg_years - focal_rev_year
    ))]

    # Protected-area screen, as in match_buffer: PA layer covering the end of the
    # BACI panel. Layers are cumulative, so PA_{focal_rev_year + 6} flags any cell
    # protected through the panel end. Applied to treated and donor plots alike.
    PA_names <- names(dplyr::select(matching_df, starts_with("PA_")))
    PA_years <- as.numeric(str_split_i(PA_names, "_", 2))
    PA_lyr_j <- PA_names[match(focal_rev_year + 6, PA_years)] %||% "PA_2016"

    matching_df_woody_veg <- matching_df %>%
      filter(!!sym(woody_veg_lyr_i) > 0 & !!sym(PA_lyr_j) < 0.5) %>%
      dplyr::select(plot_id) %>%
      collect()

    grid_points_loaded <- grid_points %>%
      filter(plot_id %in% matching_df_woody_veg$plot_id) %>%
      st_as_sf()

    grid_points_i <- inner_join(
      grid_points_loaded,
      matching_df_woody_veg,
      by = "plot_id"
    )
    grid_points_within_sf <- st_intersects_subset(
      grid_points_i,
      state_forest_mask,
      within = TRUE
    )
    focal_grid_points <- st_intersects_subset(
      grid_points_within_sf,
      focal_sf,
      negate = TRUE
    )
    focal_sf_buffer <- st_buffer(focal_sf, buffer_area_size)
    focal_grid_points <- st_intersects_subset(
      focal_grid_points,
      focal_sf_buffer,
      within = TRUE
    )
    state_forests_current <- filter(
      state_forests_clean,
      is.na(revocation_year) | revocation_year > (focal_rev_year + 5)
    )
    focal_treated <- st_difference(state_forest_mask, focal_sf) %>%
      st_union() %>%
      st_intersection(st_union(focal_sf_buffer)) %>%
      st_intersection(st_union(st_cast(
        state_forests_current,
        "MULTIPOLYGON"
      ))) %>%
      st_collection_extract("POLYGON") %>%
      st_union() %>%
      st_sf() %>%
      st_buffer(0)

    if (nrow(focal_grid_points) == 0) {
      return()
    }

    focal_grid_points <- st_intersects_subset(
      focal_grid_points,
      state_forests_current,
      within = TRUE
    )

    if (write_gpkg) {
      st_write(
        focal_sf,
        glue(
          "{OUTPUT_DIR}/focal_state_forest_gpkg/focal_state_forest_neighbour_{i}_{focal_sf_name}.gpkg"
        ),
        layer = glue("focal_state_forest_{i}_{focal_sf_name}"),
        delete_layer = TRUE,
        append = TRUE
      )
    }

    # Donor pool spatial filters for the neighbour channel:
    #   1. Must be outside the 10 km buffer of any SF revoked within the past 5 years.
    #   2. Must be inside a currently-intact state forest (not revoked within +4 years).
    #   3. Must be outside the focal SF's own 10 km buffer (not a treated unit).
    donor_pool_points <- grid_points_within_sf %>%
      st_intersects_subset(
        state_forests_clean %>%
          filter(between(
            revocation_year,
            focal_rev_year - 5,
            focal_rev_year
          )) %>%
          st_buffer(buffer_area_size),
        negate = TRUE
      ) %>%
      st_intersects_subset(
        filter(
          state_forests_clean,
          # intact through the end of the BACI panel (focal_rev_year + 5),
          # matching state_forests_current used for focal_treated above
          is.na(revocation_year) | revocation_year > (focal_rev_year + 5)
        ),
        within = TRUE
      ) %>%
      st_intersects_subset(st_buffer(focal_sf, buffer_area_size), negate = TRUE)

    if (nrow(donor_pool_points) == 0) {
      return()
    }

    matching_df_subset <- matching_df %>%
      filter(
        plot_id %in% c(focal_grid_points$plot_id, donor_pool_points$plot_id)
      ) %>%
      mutate(
        WoodyVeg = !!sym(woody_veg_lyr_i),
        WoodyVegTreat = !!sym(woody_veg_lyr_j),
        WoodyVegDelta = !!sym(woody_veg_lyr_i) - !!sym(woody_veg_lyr_j),
        PA_Final = !!sym(PA_lyr_j)
      ) %>%
      collect()

    treatment <- focal_grid_points %>%
      left_join(matching_df_subset, by = "plot_id") %>%
      mutate(
        total_clearing = rowSums(across(starts_with("slats_"))),
        PreClearing = rowSums(across(any_of(paste0(
          "slats_",
          seq(focal_rev_year - 5, focal_rev_year - 1)
        ))))
      )

    control <- donor_pool_points %>%
      left_join(matching_df_subset, by = "plot_id") %>%
      mutate(
        total_clearing = rowSums(across(starts_with("slats_"))),
        PreClearing = rowSums(across(any_of(paste0(
          "slats_",
          seq(focal_rev_year - 5, focal_rev_year - 1)
        ))))
      )

    matching_data_sf <- bind_rows(
      treatment = treatment,
      control = control,
      .id = "treat"
    ) %>%
      mutate(
        treat = treat == "treatment",
        woody_veg_start = !!sym(woody_veg_lyr_i)
      ) %>%
      filter(PA_Final < 0.5) %>%
      drop_na(treat, woody_veg_start, any_of(all.vars(formula)))

    matching_data_out <- list(
      focal_sf = focal_sf,
      focal_rev_year = focal_rev_year,
      unmatch_data = matching_data_sf,
      woody_veg_lyr = woody_veg_lyr_i
    )
    write_rds(matching_data_out, unmatch_data_path, compress = "none")
  } else {
    init(FALSE)
    matching_data_out <- read_rds(unmatch_data_path)
  }

  no_collinear <- stepwise_vif_with_cor_check(
    formula,
    st_drop_geometry(matching_data_out$unmatch_data)
  )
  safe_match <- purrr::possibly(
    function(spec) {
      do.call(
        function(...) {
          matchit_param(matching_data_out, no_collinear$formula, ...)
        },
        spec
      )
    },
    otherwise = NULL
  )
  matchit_obj <- map(spec_list, safe_match, .progress = TRUE)

  matchit_obj_path <- glue(
    "{OUTPUT_DIR}/match_data_pair_neighbour/match_data_long_all_{i}_{focal_sf_name}.rds"
  )
  write_rds(
    lapply(matchit_obj, \(x) {
      dplyr::select(
        x$match_data_long,
        any_of(c(
          "plot_id",
          "treat",
          "year",
          "BA",
          "CI",
          "t",
          "forest_clearing",
          "prop_clearing",
          "subclass",
          "woody_veg_start",
          "KoalaHabitatSuitability"
        ))
      )
    }),
    matchit_obj_path
  )

  if (!evaluate_matchit) {
    return()
  }

  matchit_summary <- matchit_list_summary(matchit_obj, spec_list)
  matchit_criteria_df <- matchit_summary %>%
    dplyr::select(any_of(c("spec_id", "ParallelTrends", "SMD", "PropKeep"))) %>%
    group_by(spec_id) %>%
    mutate(spec_id = as.numeric(spec_id)) %>%
    arrange(spec_id) %>%
    summarise(across(ParallelTrends:PropKeep, max)) %>%
    rowwise() %>%
    mutate(balance_criteria = sum(c(ParallelTrends, SMD, PropKeep)))

  matchit_criteria_idx <- matchit_criteria_df[
    matchit_criteria_df$balance_criteria ==
      max(matchit_criteria_df$balance_criteria),
    "spec_id"
  ][[1]][1]

  write_csv(
    matchit_summary,
    glue(
      "{OUTPUT_DIR}/matchit_neighbour/matchit_neighbour_summary_{i}_{focal_sf_name}.csv"
    )
  )
  write_parquet(
    matchit_obj[[matchit_criteria_idx]]$match_data_long,
    glue(
      "{OUTPUT_DIR}/match_data_pair_neighbour/match_data_long_{i}_{focal_sf_name}.parquet"
    )
  )
  if (write_gpkg) {
    st_write(
      matchit_obj[[matchit_criteria_idx]]$match_data_sp,
      glue(
        "{OUTPUT_DIR}/match_data_sp_neighbour/match_data_sp_{i}_{focal_sf_name}.gpkg"
      ),
      append = FALSE
    )
  }

  print(glue("Finished iteration {i}: {focal_sf_name}"))
  gc()
}

# =============================================================================
# Aggregate individual focal-SF geopackages into focal_sf_year.gpkg
# =============================================================================
# Run this block once after all match_neighbour() PBS jobs have completed.
# It reads every focal_state_forest_neighbour_*_Y{year}.gpkg (year-pooled groups
# only) and writes a single combined geopackage used by Scripts 03 and 04.
#
# NOTE: The woody_area_start_m2 and five_year_clearing_ha columns are set to NA
# here. They must be populated by running the woody vegetation extraction step
# (woody_veg_calc.R on the HPC) before Script 04 can compute the leakage
# multiplier or the treatment effect per hectare.
# =============================================================================
build_focal_sf_year_gpkg <- function() {
  focal_gpkg_dir <- file.path(OUTPUT_DIR, "focal_state_forest_gpkg")
  out_path <- file.path(focal_gpkg_dir, "focal_sf_year.gpkg")

  gpkg_files <- list.files(
    focal_gpkg_dir,
    pattern = "^focal_state_forest_neighbour_\\d+_Y\\d{4}\\.gpkg$",
    full.names = TRUE
  )

  if (length(gpkg_files) == 0) {
    message("No focal_state_forest_neighbour_*_Y{year}.gpkg files found in ",
            focal_gpkg_dir, " — skipping focal_sf_year.gpkg build.")
    return(invisible(NULL))
  }

  focal_list <- lapply(gpkg_files, function(f) {
    yr <- as.character(stringr::str_extract(basename(f), "\\d{4}"))
    sf_obj <- st_read(f, quiet = TRUE)
    sf_obj$run_id <- paste0("buffer_", yr)
    sf_obj$woody_area_start_m2 <- NA_real_
    sf_obj$woody_khab_area_start_m2 <- NA_real_
    sf_obj$total_area_m2 <- NA_real_
    sf_obj$woody_veg_lyr <- NA_character_
    sf_obj$five_year_clearing_ha <- NA_real_
    sf_obj[, c("run_id", "woody_area_start_m2", "woody_khab_area_start_m2",
               "total_area_m2", "woody_veg_lyr", "five_year_clearing_ha")]
  })

  focal_combined <- do.call(rbind, focal_list)
  st_write(focal_combined, out_path, delete_dsn = TRUE, quiet = TRUE)
  message("Written ", nrow(focal_combined), " rows to ", out_path)
}
