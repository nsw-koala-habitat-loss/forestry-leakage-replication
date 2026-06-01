library(INLA)

#' Calculate the Average Marginal Effect (AME) using the Incidence Rate Ratio method.
#'
#' Draws posterior samples from an INLA model fitted with a negative-binomial
#' family and an offset of log(total_woody_cells). Returns a distribution of the
#' total treatment-induced change in cleared area (in hectares) across all
#' treated plots.
#'
#' @param mdl        An INLA model object with fixed effects (Intercept), BA, CI, BACI.
#' @param n_samples  Number of posterior samples to draw.
#' @param cell_size_ha  Area of one grid cell in hectares (default 25m x 25m = 0.0625 ha).
#' @return A data.frame with columns:
#'   delta_mean_n_cells  – mean per-plot change in expected cleared cells
#'   delta_n_cells       – total change in expected cleared cells across all treated plots
#'   delta_area          – total change in cleared area (ha)
#'   prop                – delta_area as a proportion of total woody vegetation in the treated sample
generate_ame_irr <- function(
  mdl,
  n_samples  = 5000,
  cell_size_ha = 25 * 25 * 0.0001
) {
  sample <- INLA::inla.posterior.sample(
    n_samples, mdl,
    selection = list("(Intercept)" = 1, "BA" = 1, "CI" = 1, "BACI" = 1)
  )

  intercept <- sapply(sample, function(x) x$latent["(Intercept):1", 1]) %>% unname() %>% as.matrix()
  ba        <- sapply(sample, function(x) x$latent["BA:1", 1])           %>% unname() %>% as.matrix()
  ci        <- sapply(sample, function(x) x$latent["CI:1", 1])           %>% unname() %>% as.matrix()
  baci      <- sapply(sample, function(x) x$latent["BACI:1", 1])         %>% unname() %>% as.matrix()

  df <- mdl$.args$data
  treated_group    <- df[df$BACI == 1 & df$total_woody_cells > 0, ]
  total_woody_cells <- treated_group$total_woody_cells

  log_trials  <- as.matrix(rep(1, n_samples)) %*% log(total_woody_cells)
  max_counts  <- as.matrix(rep(1, n_samples)) %*% total_woody_cells

  exp_baseline <- exp(
    intercept %*% rep(1, nrow(treated_group)) +
    ba        %*% treated_group$BA +
    ci        %*% treated_group$CI +
    log_trials
  )
  exp_delta <- exp(
    intercept %*% rep(1, nrow(treated_group)) +
    ba        %*% treated_group$BA +
    ci        %*% treated_group$CI +
    baci      %*% treated_group$BACI +
    log_trials
  )

  # Cap expected counts at total woody cells
  exp_baseline <- pmin(exp_baseline, max_counts)
  exp_delta    <- pmin(exp_delta,    max_counts)

  change_in_cells <- exp_delta - exp_baseline
  delta           <- apply(change_in_cells, 1, sum)
  delta_mean      <- apply(change_in_cells, 1, mean)
  delta_area      <- delta * cell_size_ha
  total_woody     <- sum(total_woody_cells) * cell_size_ha

  data.frame(
    delta_mean_n_cells = delta_mean,
    delta_n_cells      = delta,
    delta_area         = delta_area,
    prop               = delta_area / total_woody
  )
}
