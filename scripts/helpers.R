library(ggplot2)

# ── File utilities ────────────────────────────────────────────────────────────

#' Extract run identifier and year from model output filenames.
#'
#' Matches filenames of the form  baci_model_{run}_{nn}._Y{year}.rds  and
#' returns a data.frame with columns filename, run, year.
#'
#' @param file_list  Character vector of file paths or names.
#' @param extension  File extension to match (default ".rds").
#' @param pattern    Regex with two capture groups: (run identifier) and (year).
#' @return data.frame with columns filename, run, year.
extract_years_from_files <- function(
  file_list,
  extension = ".rds",
  pattern = "_([[:alnum:]]+)_\\d+\\. Y(\\d{4})"
) {
  regex_pattern <- paste0(pattern, stringr::fixed(extension), "$")
  match_results <- stringr::str_match(file_list, regex_pattern)
  has_match <- !is.na(match_results[, 1])

  matching_files <- file_list[has_match]
  extracted_identifiers <- match_results[has_match, 2]
  extracted_years <- match_results[has_match, 3]

  if (length(matching_files) == 0) {
    message("No matching files found.")
    return(data.frame(
      filename = character(),
      run = character(),
      year = character()
    ))
  }

  data.frame(
    filename = matching_files,
    run = extracted_identifiers,
    year = extracted_years,
    stringsAsFactors = FALSE
  )
}


# ── INLA model utilities ──────────────────────────────────────────────────────

#' Load all INLA model objects from a file data.frame produced by extract_years_from_files.
pull_inla_models <- function(file_df) {
  file_df %>%
    pull(filename) %>%
    map(read_rds)
}

#' Check INLA model convergence via hyperparameter diagnostics.
#'
#' Returns TRUE if the model appears converged, FALSE otherwise.
inla_is_converged <- function(model, verbose = FALSE) {
  convergence_status <- TRUE

  if (!is.null(model$mode) && model$mode$mode.status != 0) {
    if (verbose) {
      cat("Warning: INLA optimization did not converge.\n")
    }
    convergence_status <- FALSE
  }

  hyper_summary <- model$summary.hyperpar
  if (!is.null(hyper_summary) && nrow(hyper_summary) > 0) {
    for (i in seq_len(nrow(hyper_summary))) {
      if (hyper_summary[i, "sd"] > abs(hyper_summary[i, "mean"])) {
        if (verbose) {
          cat(
            "Potential issue with hyperparameter:",
            rownames(hyper_summary)[i],
            "\n"
          )
        }
        convergence_status <- FALSE
      }
    }
  }

  # If the marginal density is still high at the integration boundary (start or end),
  # INLA may have truncated the posterior — a known issue with the negative-binomial
  # size parameter when data are sparse.
  for (marginal_name in names(model$marginals.hyperpar)) {
    marginal <- model$marginals.hyperpar[[marginal_name]]
    if (!is.null(marginal)) {
      max_y <- max(marginal[, 2])
      start_y <- marginal[1, 2]
      end_y <- marginal[nrow(marginal), 2]
      if (start_y > 0.5 * max_y || end_y > 0.5 * max_y) {
        if (verbose) {
          cat("Potential boundary truncation for:", marginal_name, "\n")
        }
        convergence_status <- FALSE
      }
    }
  }

  if (model$mode$mode.status != 0) {
    convergence_status <- FALSE
  }

  return(convergence_status)
}

#' Check that all four BACI groups (Before/After x Control/Intervention) have
#' non-zero forest clearing — a model quality guard.
data_non_zero <- function(model) {
  clearing_summary <- model$.args$data %>%
    group_by(BA, CI) %>%
    summarise(
      forest_clearing_cell = sum(forest_clearing_cell),
      .groups = "drop"
    )
  all(clearing_summary$forest_clearing_cell > 0)
}

#' Extract fixed-effect estimates and convergence diagnostics from INLA models.
#'
#' Returns a tidy data.frame with one row per model, including posterior
#' quantiles for the BACI coefficient, significance flags, sample sizes,
#' and convergence indicators.
pull_inla_estimates <- function(file_df) {
  df <- file_df %>%
    pull(filename) %>%
    map(
      function(mdl_path) {
        x <- read_rds(mdl_path)
        df <- summary(x)$fixed["BACI", ] %>% as.list()
        df$`0.01quant` <- inla.qmarginal(0.01, x$marginals.fixed$BACI)
        df$`0.99quant` <- inla.qmarginal(0.99, x$marginals.fixed$BACI)
        df$`0.05quant` <- inla.qmarginal(0.05, x$marginals.fixed$BACI)
        df$`0.95quant` <- inla.qmarginal(0.95, x$marginals.fixed$BACI)
        df$sig_pos <- df$`0.975quant` > 0 &
          df$`0.025quant` > 0 &
          df$`0.5quant` > 0
        df$sig_neg <- df$`0.975quant` < 0 &
          df$`0.025quant` < 0 &
          df$`0.5quant` < 0
        df$N <- length(unique(x$.args$data$plot_id))
        df$N_treated <- sum(x$.args$data$BACI == 1)
        df$converge <- inla_is_converged(x)
        df$all_nonzero <- data_non_zero(x)
        df
      },
      .progress = "Pulling INLA estimates"
    )

  names(df) <- paste0(pull(file_df, run), "_", pull(file_df, year))
  bind_rows(df, .id = "run_id")
}

#' Pull the focal state forest polygon from each unmatched-data RDS file.
pull_focal_sf <- function(file_df) {
  df <- file_df %>%
    pull(filename) %>%
    map(
      function(mdl_path) {
        read_rds(mdl_path)$focal_sf
      },
      .progress = "Pulling focal state forest files"
    )
  names(df) <- paste0(pull(file_df, run), "_", pull(file_df, year))
  bind_rows(df, .id = "run_id")
}


# ── MatchIt utilities ─────────────────────────────────────────────────────────

#' Load and combine MatchIt covariate-balance CSV files.
pull_matchit_summary <- function(file_df) {
  file_df$run_id <- as.character(seq_len(nrow(file_df)))
  df <- file_df %>%
    pull(filename) %>%
    map(
      read_csv,
      show_col_types = FALSE,
      progress = FALSE,
      .progress = "Pulling MatchIt Summary"
    )
  names(df) <- paste0(pull(file_df, run), "_", pull(file_df, year))
  bind_rows(df, .id = "run_id") %>%
    left_join(dplyr::select(file_df, c("run_id", "run", "year")))
}

#' Load and combine cumulative clearing CSV files.
pull_cum_clearing <- function(file_df) {
  file_df$run_id <- paste0(file_df$run, "_", file_df$year)
  df <- file_df %>%
    pull(filename) %>%
    map(
      read_csv,
      show_col_types = FALSE,
      progress = FALSE,
      .progress = "Pulling Cumulative Clearing Summary"
    )
  names(df) <- paste0(pull(file_df, run), "_", pull(file_df, year))
  bind_rows(df, .id = "run_id") %>%
    left_join(dplyr::select(file_df, c("run_id", "run", "year")))
}

#' Pull AME estimates for the best-selected matching specification.
pull_ame_all_matches <- function(file_df, matchit_summary_df) {
  df <- file_df %>%
    pull(filename) %>%
    map(
      read_csv,
      show_col_types = FALSE,
      progress = FALSE,
      .progress = "Pulling AME"
    )

  matchit_summary_df_grouped <- matchit_summary_df %>%
    filter(model_selected) %>%
    group_by(run_id) %>%
    summarise(spec_id = list(unique(spec_id)))

  name_list <- paste0(pull(file_df, run), "_", pull(file_df, year))
  names(df) <- name_list

  purrr::imap(df, \(x, idx) {
    extract_spec_id <- matchit_summary_df_grouped[
      matchit_summary_df_grouped$run_id == idx,
      "spec_id"
    ][[1]][[1]]
    unlist(x[, extract_spec_id]) %>% as.vector()
  })
}


# ── Match data utilities ──────────────────────────────────────────────────────

#' Summarise BACI clearing counts (BC, AC, BI, AI) and derived quantities
#' across the full matched dataset.
summarise_match_data <- function(df, summary_vars = c()) {
  df %>%
    mutate(
      forest_clearing_cell = ceiling(forest_clearing / 625),
      total_woody_cells = ceiling(woody_veg_start / 625),
      BACI = BA * CI
    ) %>%
    group_by(across(all_of(c(summary_vars, "plot_id", "BA", "CI")))) %>%
    summarise(
      forest_clearing_cell = sum(forest_clearing_cell),
      total_woody_cells = max(total_woody_cells),
      .groups = "drop"
    ) %>%
    mutate(
      forest_clearing_cell = if_else(
        total_woody_cells < forest_clearing_cell,
        total_woody_cells,
        forest_clearing_cell
      )
    ) %>%
    ungroup() %>%
    group_by(across(all_of(summary_vars))) %>%
    summarise(
      N = n_distinct(
        if_else(CI == 1, as.character(plot_id), NA_character_),
        na.rm = TRUE
      ),
      # BACI 2×2 cell counts: B=Before, A=After, C=Control, I=Impact (treated)
      BC = sum(forest_clearing_cell * (BA == 0) * (CI == 0)),
      AC = sum(forest_clearing_cell * (BA == 1) * (CI == 0)),
      BI = sum(forest_clearing_cell * (BA == 0) * (CI == 1)),
      AI = sum(forest_clearing_cell * (BA == 1) * (CI == 1)),
      TOT_C = sum(total_woody_cells * (CI == 0)),
      TOT_I = sum(total_woody_cells * (CI == 1)),
      .groups = "drop"
    ) %>%
    mutate(
      BACI = (AI - BI) - (AC - BC),
      BACI_prop_C = BACI / TOT_C,
      BACI_prop_I = BACI / TOT_I,
      area_ha = BACI * 625 * 0.0001
    )
}

#' Summarise BACI clearing counts at the matched-pair (subclass) level.
summarise_match_data_plot_level <- function(df) {
  df %>%
    mutate(
      forest_clearing_cell = ceiling(forest_clearing / 625),
      BACI = BA * CI
    ) %>%
    group_by(plot_id) %>%
    arrange(plot_id, t) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(total_woody_cells = ceiling(woody_veg_start / 625)) %>%
    group_by(plot_id, subclass, BA, CI) %>%
    summarise(
      forest_clearing_cell = sum(forest_clearing_cell),
      total_woody_cells = max(total_woody_cells),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      forest_clearing_cell = min(forest_clearing_cell, total_woody_cells)
    ) %>%
    ungroup() %>%
    filter(!is.na(subclass)) %>%
    group_by(subclass) %>%
    summarise(
      plot_id = first(plot_id[CI == 1]),
      BC = sum(forest_clearing_cell * (BA == 0) * (CI == 0)),
      AC = sum(forest_clearing_cell * (BA == 1) * (CI == 0)),
      BI = sum(forest_clearing_cell * (BA == 0) * (CI == 1)),
      AI = sum(forest_clearing_cell * (BA == 1) * (CI == 1)),
      TOT_C = sum(total_woody_cells * (CI == 0)),
      TOT_I = sum(total_woody_cells * (CI == 1)),
      .groups = "drop"
    ) %>%
    mutate(
      BACI = (AI - BI) - (AC - BC),
      BACI_prop_C = BACI / TOT_C,
      BACI_prop_I = BACI / TOT_I,
      area_ha = BACI * 625 * 0.0001,
      subclass = as.numeric(subclass)
    )
}

#' Compute cumulative clearing trajectories for treatment and control groups.
summarise_cum_clearing <- function(df) {
  df %>%
    mutate(forest_clearing_cell = round(forest_clearing / 625)) %>%
    group_by(plot_id) %>%
    arrange(plot_id, t) %>%
    mutate(cum_clearing = cumsum(forest_clearing_cell)) %>%
    ungroup() %>%
    group_by(CI, t) %>%
    summarise(cum_clearing = sum(cum_clearing), .groups = "drop")
}

#' Load paired match data from parquet files, selecting only required columns.
pull_match_data_pair <- function(file_df, sum_func = NULL) {
  df <- file_df %>%
    pull(filename) %>%
    map(
      function(file) {
        x <- read_parquet(
          file,
          col_select = c(
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
          )
        )
        if (!is.null(sum_func)) {
          x <- sum_func(x)
        }
        x
      },
      .progress = "Pulling match data pair"
    )
  names(df) <- paste0(pull(file_df, run), "_", pull(file_df, year))
  bind_rows(df, .id = "run_id")
}

#' Load match data across all matching specifications and apply an optional
#' summary function. Filters to specifications that meet the balance criteria.
pull_match_data_pair_all <- function(
  file_df,
  matchit_summary_df,
  sum_func = NULL
) {
  df <- file_df %>%
    pull(filename) %>%
    imap(
      function(file, idx) {
        x <- read_rds(file)
        matchit_summary_idx <- matchit_summary_df[
          matchit_summary_df$spec_id == as.numeric(idx),
        ]
        names(x) <- paste0("spec_", seq_along(x))
        x <- x[matchit_summary_idx$within_set]
        if (!is.null(sum_func)) {
          x <- map(x, sum_func, .progress = glue("Mapping summary for {idx}"))
        }
        bind_rows(x, .id = "spec_id")
      },
      .progress = "Pulling match data pair (all)"
    )
  names(df) <- paste0(pull(file_df, run), "_", pull(file_df, year))
  df
}


# ── Plotting utilities ────────────────────────────────────────────────────────

#' Rename the 'run' column from internal codes to descriptive labels for plots.
rename_run <- function(df) {
  df %>%
    mutate(
      run = factor(
        run,
        c("buffer", "neighbour"),
        c("Leakage to Private Lands", "Leakage to State Forests")
      )
    )
}
