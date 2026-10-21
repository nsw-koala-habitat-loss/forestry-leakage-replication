# =============================================================================
# Script 06: Figure 5 — per-year leakage measures + their average over years
#            + a significance grid
# =============================================================================
#   A  On-site land clearing (ha) by revocation year -- context strip
#   B  Leakage per ha of new protected area, per revocation year
#   C  Leakage per ha of new protected area, averaged over years          -- pooled, vertical
#   D  Leakage multiplier (%), per revocation year (years > 10 ha avoided clearing)
#   E  Aggregate multiplier = total leakage / total avoided clearing (%)  -- pooled, vertical
#   F  Which run-years the 95% CI excludes zero  (year x destination grid)
#
# Layout  A / B C / D E / F F  (right column narrow; A and F full width). B and C use
# one shared signed-log y-axis (linear within ~+/-0.5%, log-scaled beyond) so the
# sub-1% yearly points and the +/-50% pooled CIs read on the same scale; a couple
# of extreme years clip to a triangle/arrow. C and D share a linear scale. Pooling
# is at the posterior-draw level (each year's delta_area draws scaled to the full
# leakage zone, then combined within a draw), so year-level uncertainty propagates.
# Same delta_area x zone_woody / focal_SF_woody scaling as Figure 5 in plot_year_effects.R.
#
# Point / cell fill: white (grey in E) when the 95% CI includes zero, the
# destination colour when it excludes zero.
#
# HOW TO RUN   Same as Script 05 — needs only data/processed_data_for_plotting.RData.
# OUTPUTS (plots/)   figure5_redesign.png
# =============================================================================

library(tidyverse)
library(patchwork)
library(INLA)

source("scripts/config.R")
load(file.path(DATA_DIR, "processed_data_for_plotting.RData")) # ame_mdls, area_effect_df

standard_arrow <- arrow(length = unit(0.2, "cm"))

run_lv <- c("buffer", "neighbour")
run_labs <- c(
  buffer = "Leakage to private lands",
  neighbour = "Leakage to state forests"
)
run_short <- c(buffer = "Private\nLands", neighbour = "State\nForests")
pal <- c(buffer = "#0072B2", neighbour = "#D55E00")
q <- \(x, p) unname(quantile(x, p, na.rm = TRUE))
fill_for <- \(run, sig) if_else(sig, unname(pal[as.character(run)]), "white")

col_sc <- scale_color_manual(
  "Leakage to",
  values = pal,
  labels = c(buffer = "Private lands", neighbour = "State forests"),
  guide = guide_legend(
    override.aes = list(linetype = 0, shape = 21, fill = unname(pal))
  )
)
base <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.5),
    plot.title = element_text(size = 10, face = "bold")
  )

# ── Draw-level total leakage effect: run x year x draw, scaled to the zone ──
scaling <- area_effect_df %>%
  transmute(
    run,
    year = as.integer(year),
    ratio = woody_veg_leakage_zone_ha / woody_veg_area_ha,
    woody_veg_area_ha, # woody veg of the revoked SF = "new protected area"
    five_year_clearing_ha,
    meets_criteria = coalesce(meets_criteria, FALSE) # passed all 3 matching-QC tests
  )

tte_draws <- imap(ame_mdls, \(d, id) {
  p <- str_split_fixed(id, "_", 2)
  tibble(
    run = p[, 1],
    year = as.integer(p[, 2]),
    draw = seq_len(nrow(d)),
    tte = d$delta_area
  )
}) %>%
  bind_rows() %>%
  inner_join(scaling, by = c("run", "year")) %>%
  mutate(tte = tte * ratio)

# ── Per-year: leakage per ha of new protected area, and the multiplier ────
per_year <- tte_draws %>%
  group_by(
    run,
    year,
    five_year_clearing_ha,
    woody_veg_area_ha,
    meets_criteria
  ) %>%
  summarise(
    lb = q(tte, .025),
    med = q(tte, .5),
    ub = q(tte, .975),
    .groups = "drop"
  ) %>%
  mutate(
    run = factor(run, run_lv),
    sig = sign(lb) == sign(ub),
    pt_fill = fill_for(run, sig),
    x = year + if_else(run == "buffer", -0.22, 0.22),
    # multiplier = leakage / on-site avoided clearing
    m_lb = lb / five_year_clearing_ha,
    m_med = med / five_year_clearing_ha,
    m_ub = ub / five_year_clearing_ha
  )

# ── Pooled over years (draw level, as leakage_estimates.R::pool_effect_over_years):
#    agg_mult = aggregate multiplier = sum(leakage) / sum(avoided clearing)
pooled <- tte_draws %>%
  filter(meets_criteria) %>% # pool only years that pass all 3 matching-QC criteria
  group_by(run, draw) %>%
  summarise(
    agg_mult = sum(tte) / sum(five_year_clearing_ha),
    .groups = "drop"
  ) %>%
  group_by(run) %>%
  summarise(
    across(
      c(agg_mult),
      list(lb = \(x) q(x, .025), med = \(x) q(x, .5), ub = \(x) q(x, .975)),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) %>%
  mutate(run = factor(run, run_lv))

# ── y-scales ───────────────────────────────────────────────────────────
b_lim <- c(-2, 4) # C <-> D (shared, linear)

# ── vert_panel: a pooled estimate drawn as a vertical interval ───────────
# show_y = TRUE keeps this panel's own y-axis; y_trans/y_breaks match the
# per-year panel to its left when they share a (possibly transformed) scale.
vert_panel <- function(
  df,
  ylim,
  xlab,
  fmt = \(m) scales::comma(round(m)),
  show_y = FALSE,
  y_labels = waiver(),
  y_trans = "identity",
  y_breaks = waiver()
) {
  df <- mutate(df, x = as.integer(factor(run, run_lv)))
  ggplot(df, aes(run, med, color = run)) +
    geom_hline(
      yintercept = 0,
      linetype = 2,
      color = "gray50",
      linewidth = 0.3
    ) +
    geom_linerange(
      aes(ymin = pmax(lb, ylim[1]), ymax = pmin(ub, ylim[2])),
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    geom_segment(
      data = filter(df, ub > ylim[2]),
      aes(x = x, xend = x, y = pmax(lb, ylim[1]), yend = ylim[2]),
      arrow = standard_arrow,
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    geom_segment(
      data = filter(df, lb < ylim[1]),
      aes(x = x, xend = x, y = pmin(ub, ylim[2]), yend = ylim[1]),
      arrow = standard_arrow,
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = pt_fill),
      shape = 21,
      size = 2.7,
      stroke = 0.9,
      show.legend = FALSE
    ) +
    geom_text(
      aes(y = pmin(ub, ylim[2]), label = fmt(med)),
      vjust = -0.7,
      size = 3,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    col_sc +
    coord_cartesian(ylim = ylim, clip = "off") +
    scale_x_discrete(xlab, labels = run_short, expand = expansion(add = 0.7)) +
    scale_y_continuous(
      NULL,
      transform = y_trans,
      breaks = y_breaks,
      labels = y_labels
    ) +
    base +
    (if (show_y) {
      NULL
    } else {
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
    })
}
pct0 <- \(m) paste0(round(m * 100), "%")
pct2 <- \(m) scales::label_percent(accuracy = 0.01)(m)

# ── Context strip: on-site land clearing (ha) by revocation year ────────
# Same quantity as the multiplier denominator in C/D; x-axis aligned to the
# per-year panels so the reader can size each year's leakage against it.
clearing_df <- area_effect_df %>%
  filter(run == "buffer") %>%
  transmute(
    year = as.integer(year),
    protected_ha = woody_veg_area_ha, # woody veg added to protected areas
    avoided_ha = five_year_clearing_ha # 5-yr on-site clearing before revocation
  )
# scale factor so the avoided-clearing line shares the bar axis; undone by sec_axis
a_sf <- max(clearing_df$protected_ha, na.rm = TRUE) /
  max(clearing_df$avoided_ha, na.rm = TRUE)
avoided_col <- "#009E73"

panel_context <- ggplot(clearing_df, aes(year)) +
  geom_col(
    aes(y = protected_ha),
    fill = "gray80",
    color = "gray20",
    width = 0.65
  ) +
  geom_point(
    aes(y = avoided_ha * a_sf),
    color = 'gray90',
    size = 2,
    shape = 21,
    fill = avoided_col
  ) +
  coord_cartesian(xlim = c(1995, 2014)) +
  scale_x_continuous(NULL, breaks = seq(1995, 2014, 3)) +
  scale_y_continuous(
    "Woody veg added to\nprotected areas (ha)",
    labels = scales::label_number(scale = 1e-3, suffix = "K"),
    expand = expansion(mult = c(0, 0.1)),
    sec.axis = sec_axis(
      ~ . / a_sf,
      name = "Avoided clearing (ha)",
      labels = scales::label_comma()
    )
  ) +
  labs(title = "Protected area added and avoided clearing") +
  base +
  theme(
    axis.title.y.right = element_text(color = avoided_col),
    axis.text.y.right = element_text(color = avoided_col),
    axis.ticks.y.right = element_line(color = avoided_col),
    axis.line.y.right = element_line(color = avoided_col)
  )

# ── A / B: total leakage effect in hectares (absolute) — per year, and pooled ─
# Replaces the earlier "leakage per ha of new protected area" panels: the
# denominator normalisation was dropped in favour of showing the modelled
# effect in its own units (ha), which is what the reader ultimately cares
# about and avoids conflating "how much leaked" with "how large the new
# reserve was". Same signed-log axis trick as before (linear near zero, log
# beyond) because a couple of years/pools have CIs orders of magnitude wider
# than the rest.
a_trans <- scales::pseudo_log_trans(sigma = 50, base = 10)
a_brk <- c(-5000, -1000, -200, -50, 0, 50, 200, 1000, 5000)
a_lab <- scales::label_comma()
tr_a <- a_trans$transform
inv_a <- a_trans$inverse

# pooled sum: total leakage (ha) summed over the QC-passing years, at the
# posterior-draw level (as leakage_estimates.R::pool_effect_over_years's sum_tte)
pooled_sum <- tte_draws %>%
  filter(meets_criteria) %>%
  group_by(run, draw) %>%
  summarise(sum_tte = sum(tte), .groups = "drop") %>%
  group_by(run) %>%
  summarise(
    lb = q(sum_tte, .025),
    med = q(sum_tte, .5),
    ub = q(sum_tte, .975),
    .groups = "drop"
  ) %>%
  mutate(
    run = factor(run, run_lv),
    sig = sign(lb) == sign(ub),
    pt_fill = fill_for(run, sig)
  )

a_lim <- max(abs(c(pooled_sum$lb, pooled_sum$ub)), na.rm = TRUE) * 1.5
stub_a <- 0.28 * tr_a(a_lim)
pp <- per_year %>%
  transmute(
    run,
    x,
    pt_fill,
    lb,
    med,
    ub,
    in_rng = between(med, -a_lim, a_lim)
  ) %>%
  mutate(
    med_c = pmax(pmin(med, a_lim), -a_lim), # median clamped to the axis
    lo_c = case_when(
      lb >= a_lim ~ inv_a(tr_a(a_lim) - stub_a),
      TRUE ~ pmax(lb, -a_lim)
    ),
    hi_c = case_when(
      ub <= -a_lim ~ inv_a(-tr_a(a_lim) + stub_a),
      TRUE ~ pmin(ub, a_lim)
    )
  )

# Whiskers are always drawn (clamped to the axis); a CI tail past the axis gets
# an arrowhead, and an off-scale median is shown as a text label instead of a point.
panel_a <- ggplot(pp, aes(x, med_c, color = run)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50", linewidth = 0.3) +
  geom_linerange(
    aes(ymin = lo_c, ymax = hi_c),
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_segment(
    data = filter(pp, hi_c > a_lim | med > a_lim),
    aes(x = x, xend = x, y = lo_c, yend = a_lim),
    arrow = standard_arrow,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_segment(
    data = filter(pp, lo_c < -a_lim | med < -a_lim),
    aes(x = x, xend = x, y = hi_c, yend = -a_lim),
    arrow = standard_arrow,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_point(
    data = filter(pp, in_rng),
    aes(fill = pt_fill),
    shape = 21,
    size = 2.1,
    stroke = 0.8
  ) +
  geom_text(
    data = filter(pp, !in_rng),
    aes(
      y = if_else(med > a_lim, a_lim, -a_lim),
      label = scales::label_comma()(round(med)),
      vjust = if_else(med > a_lim, 1.8, -1.2)
    ),
    size = 2.4,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  col_sc +
  base +
  coord_cartesian(ylim = c(-a_lim, a_lim), xlim = c(1995, 2014)) +
  scale_x_continuous(NULL, breaks = seq(1995, 2014, 3)) +
  scale_y_continuous(
    "Leakage (ha)",
    transform = a_trans,
    breaks = a_brk,
    labels = a_lab
  ) +
  labs(title = "Leakage (ha)")

# B shares A's signed-log scale (same transform / breaks / limits); shows the
# QC-years summed leakage (ha), i.e. the cumulative total over the pooled years.
panel_b <- vert_panel(
  pooled_sum,
  c(-a_lim, a_lim),
  NULL,
  fmt = \(m) scales::comma(round(m)),
  show_y = TRUE,
  y_trans = a_trans,
  y_breaks = a_brk,
  y_labels = a_lab
) +
  labs(title = "Leakage sum (ha)")

# ── C: per-year leakage multiplier ─────────────────────────────────────
pb <- per_year %>%
  filter(five_year_clearing_ha > 10) %>%
  mutate(mlo = pmax(m_lb, b_lim[1]), mhi = pmin(m_ub, b_lim[2]))
panel_c <- ggplot(pb, aes(x, m_med, color = run)) +
  annotate(
    "rect",
    xmin = -Inf,
    xmax = Inf,
    ymin = -1,
    ymax = 1,
    fill = "gray95"
  ) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50", linewidth = 0.3) +
  geom_linerange(
    aes(ymin = mlo, ymax = mhi),
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  geom_point(aes(fill = pt_fill), shape = 21, size = 2.1, stroke = 0.8) +
  scale_fill_identity() +
  col_sc +
  base +
  coord_cartesian(ylim = b_lim, xlim = c(1995, 2014)) +
  scale_x_continuous(
    "Year of state forest revocation",
    breaks = seq(1995, 2014, 3)
  ) +
  scale_y_continuous(
    "Leakage multiplier (%)",
    labels = scales::label_percent()
  ) +
  labs(title = "Leakage (%)")

# ── D: aggregate multiplier (shares C's scale) ─────────────────────────
pooled_d <- pooled %>%
  transmute(
    run,
    lb = agg_mult_lb,
    med = agg_mult_med,
    ub = agg_mult_ub,
    sig = sign(lb) == sign(ub),
    pt_fill = fill_for(run, sig)
  )
panel_d <- vert_panel(
  pooled_d,
  b_lim,
  NULL,
  fmt = pct0,
  show_y = TRUE,
  y_labels = scales::label_percent()
) +
  labs(title = "Leakage (mean %)")

# ── E: significance grid (run-year x whether the 95% CI excludes zero) ──
# white cells; destination-colour fill where the 95% CI excludes zero; run-years
# that did NOT pass all 3 matching-QC criteria are diagonally hatched.
pe <- per_year %>%
  mutate(
    row = factor(
      run_labs[as.character(run)],
      levels = c(run_labs[["neighbour"]], run_labs[["buffer"]])
    ),
    cell = if_else(sig, unname(pal[as.character(run)]), "white"),
    mark = if_else(med >= 0, "+", "-")
  )

# diagonal hatch lines for the non-QC cells (built by hand — no extra package):
# 45-deg lines y = x + k, clipped to each cell's [-0.5, 0.5]^2 box.
hatch <- pe %>%
  filter(!meets_criteria) %>%
  transmute(xi = as.integer(factor(year)), yi = as.integer(row)) %>%
  distinct() %>%
  crossing(k = seq(-0.75, 0.75, by = 0.25)) %>%
  mutate(u1 = pmax(-0.5, -0.5 - k), u2 = pmin(0.5, 0.5 - k)) %>%
  filter(u1 < u2) %>%
  transmute(x = xi + u1, xend = xi + u2, y = yi + u1 + k, yend = yi + u2 + k)

panel_e <- ggplot(pe, aes(factor(year), row)) +
  geom_tile(aes(fill = cell), color = "grey80", linewidth = 0.4) +
  geom_segment(
    data = hatch,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    color = "grey55",
    linewidth = 0.3
  ) +
  geom_text(
    data = filter(pe, sig),
    aes(label = mark),
    color = "white",
    fontface = "bold",
    size = 3.8
  ) +
  geom_text(
    data = filter(pe, !sig),
    aes(label = mark),
    color = "grey45",
    size = 3
  ) +
  scale_fill_identity() +
  scale_x_discrete(NULL) +
  scale_y_discrete(NULL, labels = \(l) {
    str_to_sentence(str_remove(l, "Leakage to "))
  }) +
  labs(
    title = "Significance by revocation year"
  ) +
  base +
  theme(
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.subtitle = element_text(size = 9)
  )

# ── Compose:  A / B C / D E / F F  ───────────────────────────────────
figure5 <- panel_context +
  guide_area() +
  panel_a +
  panel_b +
  panel_c +
  panel_d +
  panel_e +
  plot_layout(
    design = "AB\nCD\nEF\nGG",
    widths = c(3, 1),
    heights = c(1, 1, 1, .5),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
  ) &
  theme(
    plot.tag.location = "panel",
    plot.tag.position = c(0.012, 0.985),
    plot.tag = element_text(hjust = 0, vjust = 1, face = "bold", size = 14),
    legend.position = "right",
    legend.direction = "vertical",
    plot.caption = element_text(hjust = 0)
  )

# F (heatmap): tag in the plot margin instead, so it doesn't sit on the tiles
figure5[[7]] <- figure5[[7]] +
  theme(plot.tag.location = "plot", plot.tag.position = "topleft")

ggsave(
  file.path(PLOT_DIR, "figure5_redesign.png"),
  figure5,
  width = 2400,
  height = 3100,
  units = "px",
  dpi = 300,
  scale = .9
)

message("Figure 5 (integrated) written to plots/figure5_redesign.png")
print(pooled)
