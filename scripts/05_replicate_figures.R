# =============================================================================
# Script 05: Replicate all paper figures
# =============================================================================
# PURPOSE
#   Generates all figures and summary tables reported in the paper.
#   This is the primary script for replication — it only requires the
#   pre-computed data file data/processed_data_for_plotting.RData.
#
# HOW TO RUN
#   1. Open an R session and set the working directory to forestry_leakage_rep/
#      (e.g. setwd("path/to/forestry_leakage_rep"))
#   2. Run this script: source("scripts/05_replicate_figures.R")
#   3. Output figures are saved to plots/
#
# INPUT
#   data/processed_data_for_plotting.RData   (pre-computed; produced by Script 04)
#   ext_data/*.xlsx                          (optional; ABARES log production data)
#
# OUTPUTS (all saved to plots/)
#   total_revocations_plt.png   — Figure 1: State forest revocations and log harvests
#   comb_coef_plts.png          — Figure 2: Main results (4-panel)
#   cum_clearing_year_plt.png   — Figure 3: Cumulative clearing by revocation year
#   focal_sf_plt.png            — Figure 4: Spatial map of leakage effects
#   comb_coef_sf_plt.png        — Figure 5: Combined dashboard
#   cumulative_avg_dist_plt.png — Figure 6: Cumulative and average effects
# =============================================================================

library(ggplot2)
library(tidyverse)
library(glue)
library(patchwork)
library(sf)
library(ggtext)
library(igraph)

has_ggmagnify <- requireNamespace("ggmagnify", quietly = TRUE)
if (has_ggmagnify) library(ggmagnify)

source("scripts/config.R")
load(file.path(DATA_DIR, "processed_data_for_plotting.RData"))

# ── Theme helper ──────────────────────────────────────────────────────────────
plot_tag_theme <- function() {
  list(theme(
    panel.background  = element_rect(color = NA, fill = "white"),
    plot.tag.position = c(0.04, 0.96),
    plot.tag.location = "panel",
    panel.border      = element_rect(color = "black", linewidth = 0.5, fill = NA)
  ))
}

scale <- 2e-1

# ── Figure 1: State forest revocations and log harvests ──────────────────────
total_revocations <- focal_sf %>%
  group_by(run_id) %>%
  summarise(
    woody_area_start_m2    = first(woody_area_start_m2),
    five_year_clearing_ha  = first(five_year_clearing_ha)
  ) %>%
  separate(run_id, c("run", "year")) %>%
  filter(run == "buffer") %>%
  mutate(total_area_ha = as.numeric(units::set_units(st_area(.), "ha"))) %>%
  st_drop_geometry() %>%
  mutate(woody_veg_start_ha = woody_area_start_m2 * 0.0001) %>%
  left_join(revocation_summary %>% mutate(year = as.character(revocation_year)), by = "year") %>%
  mutate(year = as.numeric(year)) %>%
  complete(year = full_seq(year, period = 1)) %>%
  replace_na(list(n_revocations = 0, woody_veg_start_ha = 0, total_area_ha = 0)) %>%
  left_join(abares_prod, by = join_by("year" == "Year")) %>%
  mutate(se = 0)

total_revocations_plt <- total_revocations %>%
  ggplot(aes(x = year)) +
  geom_hline(yintercept = 0, colour = "gray50") +
  geom_vline(xintercept = 1998, linetype = 2) +
  geom_vline(xintercept = 2005, linetype = 2) +
  geom_bar(aes(y = total_area_ha), stat = "identity", fill = "#e69F00") +
  geom_line(aes(y = total_prod * scale, colour = run)) +
  geom_ribbon(aes(ymin = (total_prod - se) * scale, ymax = (total_prod + se) * scale, fill = run), alpha = 0.2) +
  scale_color_viridis_d("Wood type", na.translate = FALSE) +
  scale_fill_viridis_d(guide = "none") +
  theme_bw() +
  scale_y_continuous(
    "Woody vegetation in\nrevoked area (ha)",
    labels = scales::label_number(scale = 1e-3, suffix = "K"),
    expand = expansion(mult = c(0.02, 0.1)),
    sec.axis = sec_axis(
      trans = ~ . / scale,
      name  = expression("Total Volume of Logs Harvested (m^3)"),
      labels = scales::label_number(scale = 1e-6, suffix = "M")
    )
  ) +
  coord_cartesian(xlim = c(min(total_revocations$year), max(total_revocations$year))) +
  guides(fill = "none") +
  labs(title = "State Forest Revocations and Log Harvests", x = "Year") +
  theme(panel.grid = element_blank()) +
  plot_tag_theme()

ggsave(
  "plots/total_revocations_plt.png", total_revocations_plt,
  width = 2000, height = 1500, dpi = 300, units = "px"
)

# ── Figure 2 panels: BACI coefficient and treatment effects ──────────────────

# Panel A: Raw BACI coefficient estimate
coef_plot <- baci_model_est %>%
  left_join(matchit_criteria) %>%
  rename_run() %>%
  mutate(year = as.numeric(year)) %>%
  ggplot(aes(y = `0.5quant`, x = year)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_segment(aes(y = `0.025quant`, yend = `0.975quant`, group = run, color = run),
               position = position_dodge(width = 0.5)) +
  geom_point(aes(color = run, shape = run, fill = meets_criteria),
             position = position_dodge(width = 0.5), size = 2) +
  scale_y_continuous("Difference-in-difference estimate") +
  scale_x_continuous("Year", breaks = seq(1995, 2014, 2)) +
  scale_shape_manual("", values = c(21, 22)) +
  scale_color_manual("", values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("gray70", "white")) +
  guides(fill = "none") +
  labs(x = "", title = "Model Coefficient Estimates") +
  plot_tag_theme()

if (has_ggmagnify) {
  coef_plot <- coef_plot + ggmagnify::geom_magnify(
    from = c(xmin = 1996.5, xmax = 2003, ymin = -2, ymax = 3),
    to   = c(xmin = 1997,   xmax = 2005, ymin = 15, ymax = 90),
    axes = "xy"
  )
}

# Panel B: Total treatment effect in hectares
area_effect_plt <- area_effect_df %>%
  rename_run() %>%
  mutate(year = as.numeric(year)) %>%
  ggplot(aes(y = `0.5quant`, x = year)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_segment(aes(y = `0.025quant`, yend = `0.975quant`, group = run, color = run),
               position = position_dodge(width = 0.7)) +
  geom_point(aes(color = run, fill = meets_criteria, shape = run),
             position = position_dodge(width = 0.7), size = 2) +
  scale_y_continuous("Total Treatment Effect (ha)") +
  scale_shape_manual("", values = c(21, 22)) +
  scale_color_manual("", values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("gray70", "white")) +
  coord_cartesian(xlim = c(min(total_revocations$year), max(total_revocations$year)), ylim = c(-8000, 8000)) +
  guides(fill = "none") +
  scale_x_continuous("Year", breaks = seq(1995, 2014, 2)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Treatment Effect on the Treated", x = "Year") +
  plot_tag_theme()

# Panel C: Leakage per hectare of new protected area (TTEPH)
tteph_plt <- leakage_multiplier_df %>%
  ggplot(aes(y = `tteph_0.5quant`, x = year)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_segment(aes(y = `tteph_0.025quant`, yend = `tteph_0.975quant`, group = run, color = run),
               position = position_dodge(width = 0.7)) +
  geom_point(aes(color = run, fill = meets_criteria, shape = run),
             position = position_dodge(width = 0.7), size = 2) +
  scale_x_continuous("Year", breaks = seq(1995, 2014, 2)) +
  scale_y_continuous("Leakage (per ha)") +
  scale_shape_manual("", values = c(21, 22)) +
  scale_color_manual("", values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("gray70", "white")) +
  coord_cartesian(xlim = c(min(total_revocations$year), max(total_revocations$year))) +
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  plot_tag_theme() +
  labs(title = "Leakage per hectare of new protected area")

# Panel D: Leakage multiplier (%)
leakage_mult_plt <- leakage_multiplier_df %>%
  filter(five_year_clearing_ha > 10) %>%
  ggplot(aes(y = `leak_0.5quant`, x = year)) +
  geom_rect(xmin = 1990, xmax = 2015, ymin = -1, ymax = 1, fill = "gray90") +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_segment(aes(y = `leak_0.025quant`, yend = `leak_0.975quant`, group = run, color = run),
               position = position_dodge(width = 0.7)) +
  geom_point(aes(color = run, fill = meets_criteria, shape = run),
             position = position_dodge(width = 0.7), size = 2) +
  scale_x_continuous("Year", breaks = seq(1995, 2014, 2)) +
  scale_y_continuous("Leakage (%)", labels = scales::label_percent()) +
  scale_shape_manual("", values = c(21, 22)) +
  scale_color_manual("", values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("gray70", "white")) +
  coord_cartesian(xlim = c(min(total_revocations$year), max(total_revocations$year))) +
  guides(fill = "none") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  plot_tag_theme() +
  labs(title = "Leakage multiplier")

# Combine into 4-panel Figure 2
comb_coef_plts <- (coef_plot + area_effect_plt + leakage_mult_plt + tteph_plt) +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(
  "plots/comb_coef_plts.png", comb_coef_plts,
  width = 3500, height = 2500, units = "px", dpi = 300
)

# ── Figure 3: Cumulative clearing by revocation year ─────────────────────────
cum_clearing_plt_func <- function(specified_run, years) {
  panel_letters <- letters[seq_along(years)]
  cum_clearing_all_matches_df %>%
    filter(str_detect(run_id, specified_run)) %>%
    separate(run_id, c("run", "year")) %>%
    rename_run() %>%
    mutate(criteria_score = ParallelTrends + SMD + PropKeep) %>%
    group_by(run, year) %>%
    filter(criteria_score == max(criteria_score, na.rm = TRUE)) %>%
    mutate(max_caliper = if_else(is.na(max(caliper)), NA, max(caliper))) %>%
    filter(is.na(caliper) | caliper == max_caliper) %>%
    mutate(CI = as.factor(CI)) %>%
    group_by(run, year, CI, t) %>%
    mutate(
      lb     = quantile(cum_clearing, 0),
      median = quantile(cum_clearing, 0.5),
      ub     = quantile(cum_clearing, 1)
    ) %>%
    mutate(across(c(median, lb, ub), ~ .x * 625 * 0.0001)) %>%
    mutate(
      year_annotated = if_else(criteria_score == 3, paste0("**", year, "**"), year),
      year_annotated = paste0(panel_letters[match(year, years)], ". ", year_annotated)
    ) %>%
    mutate(group = if_else(CI == 1, "Treatment", "Control")) %>%
    ggplot(aes(x = t, y = median, fill = group)) +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = 2, color = "gray50") +
    geom_line(aes(color = group), linewidth = 1) +
    facet_wrap(vars(year_annotated), scales = "free_y", ncol = 5) +
    theme_bw() +
    scale_y_continuous("Cumulative area cleared for forestry (ha)") +
    scale_x_continuous("Year from treatment year", breaks = -5:5) +
    scale_fill_manual("Group",  values = c("#CC79A7", "#0072b2")) +
    scale_color_manual("Group", values = c("#CC79A7", "#0072b2")) +
    theme(
      panel.grid      = element_blank(),
      strip.text      = element_markdown(),
      strip.background = element_rect(fill = "gray95", color = "gray50")
    )
}

years_vec <- sort(as.numeric(unique(cum_clearing_all_matches_df$year)))

cum_clearing_plt <- (
  (cum_clearing_plt_func("buffer",    years_vec) + ggtitle("Leakage to Private Lands")) /
  (cum_clearing_plt_func("neighbour", years_vec) + ggtitle("Leakage to State Forests"))
) +
  plot_layout(guides = "collect") &
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave("plots/cum_clearing_year_plt.png", cum_clearing_plt, width = 3000, height = 3000, units = "px")

# ── Figure 4 & 5: Spatial map of leakage effects ─────────────────────────────
leakage_tidy <- leakage_multiplier_df %>%
  dplyr::select(`leak_0.025quant`, `leak_0.5quant`, `leak_0.975quant`, run, year)

leakage_sf_plot_data <- focal_sf %>%
  distinct(run_id, .keep_all = TRUE) %>%
  st_transform(st_crs(nsw_border)) %>%
  separate(run_id, c("run", "year"), "_") %>%
  mutate(year = as.numeric(year)) %>%
  dplyr::select(-run) %>%
  left_join(leakage_tidy, by = "year") %>%
  filter(!is.na(`leak_0.5quant`)) %>%
  st_cast("POLYGON") %>%
  mutate(area_ha = as.numeric(units::set_units(st_area(.), "ha")), id = row_number()) %>%
  st_centroid()

neighbours_sf <- st_join(leakage_sf_plot_data, leakage_sf_plot_data,
                         join = st_is_within_distance, dist = 10000,
                         suffix = c("_orig", "_neighbour"))
edges_sf <- st_drop_geometry(dplyr::select(neighbours_sf, id_orig, id_neighbour))
graph_sf  <- graph_from_data_frame(edges_sf, directed = FALSE)
leakage_sf_plot_data$cluster_id <- as.character(components(graph_sf)$membership[
  match(leakage_sf_plot_data$id, names(components(graph_sf)$membership))
])

leakage_grouped_sf <- leakage_sf_plot_data %>%
  group_by(year, run, cluster_id) %>%
  summarise(`leak_0.5quant` = mean(`leak_0.5quant`), area_ha = sum(area_ha), .groups = "drop") %>%
  st_centroid()

leakage_sf_map <- function(specified_run, range = c(-0.02, 0.02)) {
  leakage_grouped_sf %>%
    filter(run == specified_run) %>%
    arrange(abs(`leak_0.5quant`)) %>%
    mutate(X = st_coordinates(.)[, "X"], Y = st_coordinates(.)[, "Y"]) %>%
    st_drop_geometry() %>%
    ggplot() +
    geom_sf(data = aus_border, color = "gray50", fill = "gray90") +
    geom_sf(data = nsw_border, color = "gray50", fill = "gray90") +
    geom_point(aes(fill = `leak_0.5quant`, x = X, y = Y, size = area_ha),
               color = "gray50", shape = 21) +
    coord_sf(xlim = nsw_bbox[c("xmin", "xmax")], ylim = nsw_bbox[c("ymin", "ymax")]) +
    scale_size("State Forest\nArea (ha)", range = c(1, 4)) +
    scale_fill_gradient2("Leakage", low = "#0072b2", high = "#d55e00",
                         limits = range, oob = scales::squish, labels = scales::percent) +
    guides(size = "none") +
    theme_void() +
    plot_tag_theme() +
    theme(legend.position = "bottom", legend.key.width = unit(0.7, "null"), plot.title.position = "plot")
}

leakage_sf_buffer_plt   <- leakage_sf_map("buffer",    c(-1.1, 1.1)) + labs(title = "Leakage effect to Private Lands")
leakage_sf_neighbour_plt <- leakage_sf_map("neighbour", c(-1.1, 1.1)) + labs(title = "Leakage effect to State Forests")

leakage_sf_comb_plt <- leakage_sf_buffer_plt + leakage_sf_neighbour_plt +
  plot_annotation(tag_levels = "A") &
  plot_layout(guides = "keep", nrow = 1)

ggsave("plots/focal_sf_plt.png", leakage_sf_comb_plt, dpi = 300, width = 3200, height = 1700, units = "px")

# Combined dashboard (Figure 5)
comb_coef_sf_plt <- (
  (total_revocations_plt + area_effect_plt +
    (leakage_sf_buffer_plt + plot_layout(guides = "keep")) +
    coef_plot + leakage_mult_plt +
    (leakage_sf_neighbour_plt + plot_layout(guides = "keep"))) +
    plot_layout(widths = c(0.5, 0.5, 0.6), nrow = 2, guides = "collect") &
    plot_annotation(tag_levels = "A")
) & theme(legend.position = "bottom")

ggsave("plots/comb_coef_sf_plt.png", comb_coef_sf_plt, units = "px", height = 2700, width = 4500)

# ── Summary tables ────────────────────────────────────────────────────────────

revocation_area_summary <- woody_veg_area %>%
  dplyr::select(run, year, woody_veg_leakage_zone_ha) %>%
  left_join(focal_sf_tidy %>% mutate(year = as.numeric(year))) %>%
  left_join(sample_size %>% separate(run_id, c("run","year")) %>% mutate(year = as.numeric(year))) %>%
  rename_run() %>%
  mutate(woody_veg_focal_sf_ha = woody_area_start_m2 * 0.0001) %>%
  dplyr::select(run, year, n, woody_veg_focal_sf_ha, woody_veg_leakage_zone_ha) %>%
  mutate(across(woody_veg_focal_sf_ha:woody_veg_leakage_zone_ha, round))

match_data_pair_summarised_range <- match_data_pair_all_summarised %>%
  dplyr::select(-any_of(contains(c("prop","area")))) %>%
  group_by(run_id) %>%
  mutate(across(-any_of(c("run_id","spec_id")), ~ .x * 25 * 25 * 0.0001)) %>%
  reframe(across(-spec_id, \(x) quantile(x, p = c(0.025, 0.5, 0.975)))) %>%
  mutate(var = rep(c("lb","median","ub"), length.out = nrow(.))) %>%
  mutate(across(-any_of(c("run_id","spec_id","var")), ~ round(.x, 1))) %>%
  pivot_wider(id_cols = run_id, values_from = -all_of(c("run_id","var")), names_from = var) %>%
  separate(run_id, c("run","year")) %>%
  rename_run() %>%
  mutate(
    BC_range   = glue("[{BC_lb} - {BC_ub}]"),
    AC_range   = glue("[{AC_lb} - {AC_ub}]"),
    BACI_range = glue("[{BACI_lb} - {BACI_ub}]")
  ) %>%
  dplyr::select(run, year, BC_median, BC_range, BI_median, AC_median, AC_range, AI_median, BACI_median, BACI_range)

# ── Export tidy spatial leakage estimates ─────────────────────────────────────
focal_sf_tidy_geom <- focal_sf %>%
  mutate(woody_veg_area_ha = woody_area_start_m2 * 0.0001) %>%
  distinct(run_id, .keep_all = TRUE) %>%
  separate(run_id, c("run","year")) %>%
  dplyr::select(-run)

leakage_sf_tidy <- focal_sf_tidy_geom %>%
  left_join(area_effect_df) %>%
  mutate(year = as.numeric(year), TTE_lb = `0.025quant`, TTE_median = `0.5quant`, TTE_ub = `0.975quant`) %>%
  dplyr::select(any_of(c("run","year","woody_veg_area_ha","TTE_lb","TTE_median","TTE_ub"))) %>%
  rename_run() %>%
  left_join(dplyr::select(leakage_multiplier_df, -woody_veg_area_ha)) %>%
  mutate(
    leakage_lb      = `leak_0.025quant`,
    leakage_median  = `leak_0.5quant`,
    leakage_ub      = `leak_0.975quant`,
    TTE_perha_lb    = `tteph_0.025quant`,
    TTE_perha_median = `tteph_0.5quant`,
    TTE_perha_ub    = `tteph_0.975quant`
  ) %>%
  mutate(across(leakage_lb:leakage_ub, ~ if_else(is.infinite(.), NA, .))) %>%
  dplyr::select(any_of(c("run","year","woody_veg_area_ha",
                          "TTE_lb","TTE_median","TTE_ub",
                          "leakage_lb","leakage_median","leakage_ub",
                          "TTE_perha_lb","TTE_perha_median","TTE_perha_ub")))

leakage_out_path <- file.path(DATA_DIR, "leakage_sf_tidy.gpkg")
if (file.exists(leakage_out_path)) file.remove(leakage_out_path)
st_write(leakage_sf_tidy, leakage_out_path)

# ── Figure 6: Cumulative average and total treatment effects ─────────────────
sf_area_info <- woody_veg_area %>%
  rename_run() %>%
  left_join(leakage_sf_tidy) %>%
  dplyr::select(run, year, woody_veg_area_ha, woody_veg_leakage_zone_ha)

tte_dist <- ame_mdls %>%
  map(~ pull(., "delta_area")) %>%
  bind_rows(.id = "run") %>%
  mutate(id = row_number()) %>%
  pivot_longer(cols = -id, values_to = "sample_effect") %>%
  left_join(woody_veg_ha_sample, by = join_by("name" == "run_id")) %>%
  mutate(tte_prop = sample_effect / woody_veg_ha) %>%
  separate(name, c("run","year")) %>%
  rename_run() %>%
  mutate(year = as.numeric(year)) %>%
  left_join(dplyr::select(sf_area_info, c("woody_veg_area_ha","woody_veg_leakage_zone_ha","run","year"))) %>%
  mutate(tteph = tte_prop * woody_veg_leakage_zone_ha / woody_veg_area_ha,
         tte   = tte_prop * woody_veg_leakage_zone_ha) %>%
  group_by(run, id) %>%
  left_join(mutate(rename_run(matchit_criteria), year = as.numeric(year)))

tte_year <- tte_dist %>%
  group_by(run, id) %>%
  arrange(year) %>%
  mutate(
    tteph_cum = cumsum(tteph * woody_veg_area_ha) / cumsum(woody_veg_area_ha),
    tte_cum   = cumsum(sample_effect)
  ) %>%
  group_by(run, year) %>%
  summarise(
    tteph_cum_lb     = quantile(tteph_cum, 0.025),
    tteph_cum_median = quantile(tteph_cum, 0.5),
    tteph_cum_ub     = quantile(tteph_cum, 0.975),
    tte_cum_lb       = quantile(tte_cum, 0.025),
    tte_cum_median   = quantile(tte_cum, 0.5),
    tte_cum_ub       = quantile(tte_cum, 0.975)
  )

ylim_cum     <- c(-0.025, 0.1)
ylim_cum_tte <- c(-2200, 1300)

cumulative_plt <- ggplot(tte_year, aes(x = year, y = tteph_cum_median)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_ribbon(aes(ymin = tteph_cum_lb, ymax = tteph_cum_ub, group = run, fill = run), alpha = 0.3) +
  geom_point(aes(color = run), size = 2) +
  geom_line(aes(color = run)) +
  scale_fill_manual("",  values = c("#0072B2", "#D55E00")) +
  scale_color_manual("", values = c("#0072B2", "#D55E00")) +
  scale_y_continuous("Total Treatment\nEffect per ha") +
  scale_x_continuous("Year", breaks = seq(1995, 2014, 2)) +
  coord_cartesian(ylim = ylim_cum) +
  ggtitle("Cumulative Weighted Average Effect") +
  theme_bw() + theme(panel.grid = element_blank())

cumulative_tte_plt <- ggplot(tte_year, aes(x = year, y = tte_cum_median)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_ribbon(aes(ymin = tte_cum_lb, ymax = tte_cum_ub, group = run, fill = run), alpha = 0.3) +
  geom_point(aes(color = run), size = 2) +
  geom_line(aes(color = run)) +
  scale_fill_manual("",  values = c("#0072B2", "#D55E00")) +
  scale_color_manual("", values = c("#0072B2", "#D55E00")) +
  scale_y_continuous("Total Treatment\nEffect") +
  scale_x_continuous("Year", breaks = seq(1995, 2014, 2)) +
  coord_cartesian(ylim = ylim_cum_tte) +
  ggtitle("Cumulative Total Treatment Effect") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())

avg_dist <- tte_dist %>%
  summarise(tteph = sum(tteph * woody_veg_area_ha) / sum(woody_veg_area_ha),
            tte   = sum(sample_effect))

avg_dist_plt <- avg_dist %>%
  ggplot(aes(y = tteph)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_density(aes(color = run)) +
  geom_density(aes(fill = run), alpha = 0.2, color = NA) +
  ylim(ylim_cum) +
  scale_color_manual("", values = c("#0072B2", "#D55E00")) +
  scale_fill_manual("",  values = c("#0072B2", "#D55E00")) +
  theme_bw() + guides(color = "none", fill = "none") +
  scale_x_continuous("Density") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggtitle("Average Effect per ha")

avg_dist_tte_plt <- avg_dist %>%
  ggplot(aes(y = tte)) +
  geom_hline(yintercept = 0, color = "gray50", linetype = 2) +
  geom_density(aes(color = run)) +
  geom_density(aes(fill = run), alpha = 0.2, color = NA) +
  ylim(ylim_cum_tte) +
  scale_color_manual("", values = c("#0072B2", "#D55E00")) +
  scale_fill_manual("",  values = c("#0072B2", "#D55E00")) +
  theme_bw() + guides(color = "none", fill = "none") +
  scale_x_continuous("") +
  theme(panel.grid = element_blank(),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  ggtitle("Total Treatment Effect")

cumulative_avg_dist_plt <- (cumulative_tte_plt + avg_dist_tte_plt + cumulative_plt + avg_dist_plt) +
  plot_layout(guides = "collect", widths = c(0.7, 0.3)) &
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(
  "plots/cumulative_avg_dist_plt.png", cumulative_avg_dist_plt,
  dpi = 300, width = 2500, height = 2500, units = "px"
)

message("All figures saved to plots/")
