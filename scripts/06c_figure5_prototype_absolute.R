# =============================================================================
# Script 06c: Prototype panels for showing the leakage effect in absolute ha
#             (candidates to replace Script 06's per-year panel A)
# =============================================================================
# Two alternatives to the "total leakage effect (ha), by year" panel in
# Script 06, both keeping the effect in ha rather than normalising by any
# denominator:
#
#   Prototype 1 (cumulative)   Running total of pooled QC-year leakage (ha) as
#                               revocations accrue, drawn like a cumulative
#                               dose-response curve. Pooling is at the
#                               posterior-draw level (cumsum within each draw,
#                               then quantiles across draws), so year-to-year
#                               correlation in the uncertainty is preserved.
#   Prototype 2 (side-by-side) Per revocation year, the avoided on-site
#                               clearing (ha) directly beside the modelled
#                               leakage effect (ha) for each destination, so
#                               the reader can size one against the other
#                               without a ratio.
#
# Neither is wired into figure5_redesign.png yet -- this script renders both
# side by side in one comparison PNG so a choice can be made, matching Script
# 06's per-year draw-level pooling (tte_draws) exactly.
#
# HOW TO RUN   Same as Script 06 -- needs only data/processed_data_for_plotting.RData.
# OUTPUTS (plots/)   figure5_absolute_prototypes.png
# =============================================================================

library(tidyverse)
library(patchwork)

source("scripts/config.R")
load(file.path(DATA_DIR, "processed_data_for_plotting.RData")) # ame_mdls, area_effect_df

run_lv <- c("buffer", "neighbour")
run_labs <- c(buffer = "Private lands", neighbour = "State forests")
pal <- c(buffer = "#0072B2", neighbour = "#D55E00")
q <- \(x, p) unname(quantile(x, p, na.rm = TRUE))

base <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.5),
    plot.title = element_text(size = 10, face = "bold")
  )
col_sc <- scale_color_manual("Leakage to", values = pal, labels = run_labs)
fill_sc <- scale_fill_manual("Leakage to", values = pal, labels = run_labs)

# ── Draw-level total leakage effect: run x year x draw, scaled to the zone ──
scaling <- area_effect_df %>%
  transmute(
    run,
    year = as.integer(year),
    ratio = woody_veg_leakage_zone_ha / woody_veg_area_ha,
    five_year_clearing_ha,
    meets_criteria = coalesce(meets_criteria, FALSE)
  )

tte_draws <- imap(ame_mdls, \(d, id) {
  p <- str_split_fixed(id, "_", 2)
  tibble(run = p[, 1], year = as.integer(p[, 2]), draw = seq_len(nrow(d)), tte = d$delta_area)
}) %>%
  bind_rows() %>%
  inner_join(scaling, by = c("run", "year")) %>%
  mutate(tte = tte * ratio)

# =============================================================================
# Prototype 1: cumulative absolute leakage (ha) over years, QC years only,
# pooled at the draw level -- cumsum(tte) within each draw across years in
# chronological order, then take quantiles of the cumulative total across
# draws at each year.
# =============================================================================
cum_draws <- tte_draws %>%
  filter(meets_criteria) %>%
  arrange(run, draw, year) %>%
  group_by(run, draw) %>%
  mutate(cum_tte = cumsum(tte)) %>%
  ungroup()

cum_summary <- cum_draws %>%
  group_by(run, year) %>%
  summarise(lb = q(cum_tte, .025), med = q(cum_tte, .5), ub = q(cum_tte, .975), .groups = "drop") %>%
  mutate(run = factor(run, run_lv))

panel_cumulative <- ggplot(cum_summary, aes(year, med, color = run, fill = run)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  col_sc +
  fill_sc +
  base +
  scale_x_continuous("Revocation year", breaks = scales::pretty_breaks()) +
  scale_y_continuous("Cumulative leakage effect (ha)", labels = scales::label_comma()) +
  labs(
    title = "Prototype 1: Cumulative absolute leakage over years",
    subtitle = "Running total across QC-passing revocation years, pooled at the posterior-draw level"
  ) +
  theme(legend.position = "bottom")

# =============================================================================
# Prototype 2: avoided clearing (ha) vs. leakage effect (ha), side by side,
# per revocation year -- no ratio taken, so the reader compares magnitudes
# directly on one shared ha axis.
# =============================================================================
per_year <- tte_draws %>%
  group_by(run, year, five_year_clearing_ha, meets_criteria) %>%
  summarise(lb = q(tte, .025), med = q(tte, .5), ub = q(tte, .975), .groups = "drop") %>%
  mutate(run = factor(run, run_lv))

clearing_bar <- per_year %>%
  distinct(year, five_year_clearing_ha) %>%
  mutate(run = "Avoided clearing\n(on-site)")

panel_sidebyside <- ggplot(per_year, aes(x = year)) +
  geom_col(
    data = clearing_bar,
    aes(y = five_year_clearing_ha),
    fill = "gray85",
    color = "gray50",
    width = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50", linewidth = 0.3) +
  geom_pointrange(
    aes(y = med, ymin = lb, ymax = ub, color = run),
    position = position_dodge(width = 1.4),
    size = 0.35,
    fatten = 2.5
  ) +
  col_sc +
  base +
  scale_x_continuous("Revocation year", breaks = scales::pretty_breaks()) +
  scale_y_continuous("Hectares", labels = scales::label_comma()) +
  coord_cartesian(ylim = c(-2000, 2000)) +
  labs(
    title = "Prototype 2: Avoided clearing vs. leakage effect, both in ha",
    subtitle = "Grey bars = 5-yr avoided on-site clearing; points = modelled leakage effect (95% CI), by destination"
  ) +
  theme(legend.position = "bottom")

comb <- panel_cumulative / panel_sidebyside +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(PLOT_DIR, "figure5_absolute_prototypes.png"),
  comb,
  width = 2000,
  height = 2400,
  units = "px",
  dpi = 300
)

message("Prototype panels written to plots/figure5_absolute_prototypes.png")
