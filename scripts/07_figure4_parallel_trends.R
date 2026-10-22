# =============================================================================
# Script 07: Figure 4 — parallel-trends check (re-based)
# =============================================================================
# Added for the R&R. Needs only data/processed_data_for_plotting.RData.
# Run from the package root:  source("scripts/07_figure4_parallel_trends.R").
#
# Both treatment and control arms are set to 0 in the year before revocation
# (t = -1) and shown as cumulative clearing from that point, as a % of woody
# vegetation -- the quantity the parallel-trends test evaluates. The assumption
# holds when the two lines coincide inside the shaded pre-revocation window; a
# gap after t = 0 is the leakage effect.
#
# Shaded band = control line +/- 1.96 x matched-pair SE (SE of the within-pair
# treatment - control difference). Where the treatment line sits inside it, the
# groups are not distinguishable that year. Red panel labels: the selected
# matched specification fails the pre-trend criterion (parallel_trends p >= 0.05).
#
# Output:  plots/fig4_parallel_trends.png
# =============================================================================

library(tidyverse)
library(patchwork)
library(ggtext)

source("scripts/config.R")
load(file.path(DATA_DIR, "processed_data_for_plotting.RData"))

# ---- pass / fail of the pre-trend test for the selected matched spec ---------
pt <- matchit_summary_df %>%
  filter(model_selected) %>%
  separate(run_id, c("run", "year"), remove = FALSE) %>%
  select(run, year, pass = ParallelTrends)

n_pairs <- match_data_pair_full %>%
  distinct(run_id, plot_id, CI) %>%
  count(run_id, CI) %>%
  pivot_wider(names_from = CI, values_from = n, names_prefix = "n") %>%
  separate(run_id, c("run", "year")) %>%
  transmute(run, year, n_lab = sprintf("%s/%s", n1, n0))

# ---- re-base each plot's cumulative clearing proportion at t = -1 ------------
rel <- match_data_pair_full %>%
  arrange(run_id, plot_id, t) %>%
  group_by(run_id, plot_id) %>%
  mutate(cum = cumsum(prop_clearing)) %>%
  group_by(run_id, plot_id, CI) %>%
  mutate(rel = 100 * (cum - cum[t == -1])) %>%
  ungroup()

arm_mean <- rel %>%
  group_by(run_id, CI, t) %>%
  summarise(m = mean(rel), .groups = "drop")

# matched-pair SE: within each matched set (subclass), mean(treated) - mean(control)
pair_se <- rel %>%
  group_by(run_id, subclass, CI, t) %>%
  summarise(rel = mean(rel), .groups = "drop") %>%
  pivot_wider(names_from = CI, values_from = rel, names_prefix = "a") %>%
  filter(!is.na(a0), !is.na(a1)) %>%
  group_by(run_id, t) %>%
  summarise(se_pair = sd(a1 - a0) / sqrt(n()), .groups = "drop")

plot_df <- arm_mean %>%
  left_join(pair_se, by = c("run_id", "t")) %>%
  separate(run_id, c("run", "year")) %>%
  left_join(pt, by = c("run", "year")) %>%
  left_join(n_pairs, by = c("run", "year")) %>%
  mutate(
    group = if_else(CI == 1, "Treatment", "Control"),
    lo = if_else(group == "Control", m - 1.96 * se_pair, NA_real_),
    hi = if_else(group == "Control", m + 1.96 * se_pair, NA_real_)
  )

# ---- one channel panel ------------------------------------------------------
pt_panel <- function(channel, ttl) {
  d <- plot_df %>% filter(run == channel)
  yrs <- sort(as.integer(unique(d$year)))
  d <- d %>%
    mutate(
      idx = match(as.integer(year), yrs),
      lab = if_else(
        pass,
        sprintf("(%d) **%s**", idx, year),
        sprintf(
          "(%d) <span style='color:#c0392b'>**%s** &mdash; fails</span>",
          idx,
          year
        )
      ),
      lab = fct_reorder(lab, idx)
    )
  n_lab_df <- distinct(d, year, lab, n_lab)

  ggplot(d, aes(t, m, colour = group, fill = group)) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = -1,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey93"
    ) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    #geom_ribbon(aes(ymin = lo, ymax = hi), colour = NA, alpha = 0.3, na.rm = TRUE) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.1) +
    geom_text(
      data = n_lab_df,
      aes(x = -Inf, y = Inf, label = paste0("n ", n_lab)),
      inherit.aes = FALSE,
      hjust = -0.12,
      vjust = 1.5,
      size = 2.5,
      colour = "grey45"
    ) +
    facet_wrap(vars(lab), scales = "free_y", ncol = 3) +
    scale_colour_manual(
      "Group",
      values = c(Control = "#CC79A7", Treatment = "#0072b2")
    ) +
    scale_fill_manual(
      "Group",
      values = c(Control = "#CC79A7", Treatment = "#0072b2")
    ) +
    scale_x_continuous("Year from revocation", breaks = -5:5) +
    scale_y_continuous(
      "Cumulative clearing since the year before revocation  (% of woody vegetation)",
      labels = function(x) paste0(x, "%")
    ) +
    ggtitle(ttl) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      strip.text = element_markdown(size = 9),
      strip.background = element_rect(fill = "gray95", colour = "gray50"),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

fig4_pt <- (pt_panel("buffer", "A. Leakage to private lands") |
  pt_panel("neighbour", "B. Leakage to state forests")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig4_pt <- fig4_pt +
  plot_annotation(
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30")
    )
  )

ggsave(
  "plots/fig4_parallel_trends.png",
  fig4_pt,
  width = 3800,
  height = 2550,
  units = "px"
)
