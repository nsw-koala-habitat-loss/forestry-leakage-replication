# =============================================================================
# Script 06b: Derivation table for the leakage multiplier (+ year x channel heatmap)
# =============================================================================
# Companion to Script 06. Shows, step by step, how the per-year leakage estimate
# is built:  raw matched-sample clearing counts (BC/AC/BI/AI)  ->  diff-in-diff
# (BACI)  ->  matching-QC status  ->  model AME scaled to the full leakage zone
# (Step 3)  ->  leakage multiplier (Step 4); then the pooled-over-years estimates.
#
#   Matching QC  the three balance criteria on the model-selected matching spec
#                (parallel pre-trends p >= 0.05; SMD < 0.25; proportion of
#                treated units retained >= 0.8) and their conjunction
#                `meets_criteria`.  Figure 5 pools only years where all 3 hold.
#   Pooled       per-year effects combined at the posterior-draw level (as in
#                Script 06 and leakage_estimates.R::pool_effect_over_years): each
#                year's delta_area draws are scaled to the leakage zone (Step 3),
#                then within each of 10,000 draws summed across years (Step 3)
#                or combined as sum(leakage)/sum(avoided clearing) (Step 4
#                aggregate); quantiles are then taken over draws.  Shown for the
#                QC-passing years (the Figure 5 pooling) and for all fitted years.
#
# Uses Script 04's leakage_multiplier_df + match_data_pair_df + area_effect_df +
# ame_mdls unchanged, so the numbers match Figure 5.  The per-year Step 3 bounds
# are quantiles of delta_area times the (deterministic, positive) zone/focal-SF
# area ratio -- identical to Script 06's scale-then-quantile order.
#
# HOW TO RUN   Same as Script 05 — needs only data/processed_data_for_plotting.RData.
#             Requires {gt}; {openxlsx} (or {writexl}) for the .xlsx export.
# OUTPUTS (plots/)
#   figure5_derivation_table.html   — always written
#   figure5_derivation_table.png    — if a headless browser is available
#   figure5_derivation_table.xlsx   — formatted (openxlsx) or plain (writexl) workbook
#   figure5_derivation_heatmap.png  — leakage multiplier by year x destination
# =============================================================================

library(tidyverse)
library(scales)
library(gt)

source("scripts/config.R")
load(file.path(DATA_DIR, "processed_data_for_plotting.RData"))

run_levels <- c("Leakage to Private Lands", "Leakage to State Forests")
qc_lab <- function(x) dplyr::case_when(x ~ "pass", !x ~ "fail", TRUE ~ NA_character_)

# ── Per-year effect table, from Script 04's leakage_multiplier_df ────────────
# tte_* = total leakage effect (ha); leak_* = effect / on-site avoided clearing.
# `run` is already the descriptive factor (Script 04 calls rename_run()).
effect_df <- leakage_multiplier_df %>%
  mutate(sig = sig_pos | sig_neg) %>%
  transmute(
    run = factor(run, run_levels),
    year = as.numeric(year),
    tte_lb = `0.025quant`, tte_med = `0.5quant`, tte_ub = `0.975quant`,
    leak_lb = `leak_0.025quant`, leak_med = `leak_0.5quant`, leak_ub = `leak_0.975quant`,
    sig
  )

# Matching-QC status *and* the new-protected-area size, for every run-year that
# has a matched sample -- area_effect_df is the full 29-run set, before the
# all_nonzero filter that trims effect_df/leakage_multiplier_df. Sourcing
# woody_veg_area_ha from here (rather than from effect_df) means a run-year
# still shows its protected-area size even when all_nonzero == FALSE drops it
# from the leakage-effect columns (e.g. 2014, both channels); a genuinely
# missing area (e.g. Private Lands 2012, NA upstream in focal_sf_year.gpkg --
# no model was ever fit for that run-year) still shows "--".
qc_df <- area_effect_df %>%
  mutate(year = as.numeric(year)) %>%
  rename_run() %>%
  transmute(
    run, year,
    woody_veg_area_ha, five_year_clearing_ha,
    ParallelTrends, SMD, PropKeep,
    meets_criteria = coalesce(meets_criteria, FALSE)
  )

# ── Draw-level pooling: rebuild Script 06's tte_draws, then pool ────────────
scaling <- area_effect_df %>%
  transmute(
    run, year = as.integer(year),
    ratio = woody_veg_leakage_zone_ha / woody_veg_area_ha,
    woody_veg_area_ha, five_year_clearing_ha,
    meets_criteria = coalesce(meets_criteria, FALSE)
  )

tte_draws <- imap(ame_mdls, \(d, id) {
  p <- str_split_fixed(id, "_", 2)
  tibble(run = p[, 1], year = as.integer(p[, 2]), draw = seq_len(nrow(d)), tte = d$delta_area)
}) %>%
  bind_rows() %>%
  inner_join(scaling, by = c("run", "year")) %>%
  mutate(tte = tte * ratio)

# pool `draws` within each posterior draw, then take quantiles over draws.
pool_rows <- function(draws, tag) {
  draws %>%
    group_by(run, draw) %>%
    summarise(
      n_years    = dplyr::n_distinct(year),
      sum_tte    = sum(tte),                               # Step 3, pooled (total over years)
      agg_mult   = sum(tte) / sum(five_year_clearing_ha),  # Step 4, aggregate (sum / sum)
      .groups = "drop"
    ) %>%
    group_by(run) %>%
    summarise(
      n_years   = first(n_years),
      tte_lb    = quantile(sum_tte, 0.025),
      tte_med   = quantile(sum_tte, 0.5),
      tte_ub    = quantile(sum_tte, 0.975),
      leak_med  = quantile(agg_mult, 0.5),
      .groups = "drop"
    ) %>%
    left_join(
      draws %>%
        distinct(run, year, woody_veg_area_ha, five_year_clearing_ha) %>%
        group_by(run) %>%
        summarise(
          woody_veg_area_ha = sum(woody_veg_area_ha, na.rm = TRUE),
          five_year_clearing_ha = sum(five_year_clearing_ha, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "run"
    ) %>%
    mutate(
      run  = factor(run, c("buffer", "neighbour"), run_levels),
      year = sprintf("Pooled %s (%d)", tag, n_years),
      sig  = if_else(sign(tte_lb) == sign(tte_ub), "significant", "n.s.")
    )
}

pooled_df <- bind_rows(
  pool_rows(dplyr::filter(tte_draws, meets_criteria), "QC"),   # <- the pooling used in Figure 5
  pool_rows(tte_draws, "all")
)

# ── Derivation table: raw counts -> QC -> AME -> multiplier -> per-ha ────────
per_year_df <- match_data_pair_df %>%
  separate(run_id, c("run", "year"), sep = "_") %>%
  mutate(year = as.numeric(year)) %>%
  rename_run() %>%                       # match effect_df's renamed `run` before the join
  left_join(
    dplyr::select(effect_df, run, year, sig, tte_lb, tte_med, tte_ub),
    by = c("run", "year")
  ) %>%
  left_join(qc_df, by = c("run", "year")) %>%
  mutate(
    leak_med = tte_med / five_year_clearing_ha,   # Step 4: leakage multiplier
    BACI_ha = BACI * 625 * 0.0001,
    sig = if_else(sig, "significant", "n.s."),
    meets_criteria = coalesce(meets_criteria, FALSE),
    QC = if_else(meets_criteria, "yes", "no"),
    dplyr::across(c(ParallelTrends, SMD, PropKeep), qc_lab)
  ) %>%
  arrange(run, year) %>%
  mutate(year = as.character(year))

derivation_df <- bind_rows(per_year_df, pooled_df) %>%
  arrange(run)   # stable: per-year rows (by year) then the two pooled rows, per run

derivation_table <- derivation_df %>%
  dplyr::select(run, year, BC, AC, BI, AI, BACI, BACI_ha,
                ParallelTrends, SMD, PropKeep, QC,
                tte_lb, tte_med, tte_ub, woody_veg_area_ha,
                five_year_clearing_ha, leak_med, sig) %>%
  gt(groupname_col = "run", rowname_col = "year") %>%
  tab_header(
    title = "Derivation of the leakage multiplier",
    subtitle = "Matched-sample clearing -> total leakage effect -> multiplier, with matching-QC status and the pooled-over-years estimates"
  ) %>%
  tab_spanner("Step 1: matched-sample clearing (grid cells)", c(BC, AC, BI, AI)) %>%
  tab_spanner("Step 2: diff-in-diff", c(BACI, BACI_ha)) %>%
  tab_spanner("Matching QC (model-selected spec)", c(ParallelTrends, SMD, PropKeep, QC)) %>%
  tab_spanner("Step 3: total leakage effect (ha, scaled to full leakage zone)",
              c(tte_lb, tte_med, tte_ub, woody_veg_area_ha)) %>%
  tab_spanner("Step 4: leakage multiplier", c(five_year_clearing_ha, leak_med)) %>%
  cols_label(
    BC = "Before, control", AC = "After, control",
    BI = "Before, treated", AI = "After, treated",
    BACI = "BACI (cells)", BACI_ha = "BACI (ha)",
    ParallelTrends = "Parallel trends", SMD = "SMD", PropKeep = "Prop. kept", QC = "All 3 met",
    tte_lb = "2.5%", tte_med = "Median", tte_ub = "97.5%",
    woody_veg_area_ha = "New protected area (ha woody veg.)",
    five_year_clearing_ha = "Avoided clearing (ha)", leak_med = "Multiplier",
    sig = "95% CI"
  ) %>%
  fmt_number(c(BC, AC, BI, AI, BACI), decimals = 0) %>%
  fmt_number(c(BACI_ha, tte_lb, tte_med, tte_ub, woody_veg_area_ha, five_year_clearing_ha),
             decimals = 0) %>%
  fmt_percent(leak_med, decimals = 0) %>%
  data_color(leak_med, palette = c("#0072B2", "white", "#D55E00"), domain = c(-2, 2)) %>%
  sub_missing(columns = everything(), missing_text = "--") %>%
  tab_style(cell_text(style = "italic"),
            locations = cells_body(rows = grepl("^Pooled", year))) %>%
  tab_footnote("BACI = (After, treated - Before, treated) - (After, control - Before, control); a cell is 625 m^2 (25 x 25 m).") %>%
  tab_footnote("Matching QC: parallel pre-trends (p >= 0.05), SMD (< 0.25), proportion of treated units retained (>= 0.8), on the model-selected matching spec. 'All 3 met' = meets_criteria; '--' where no matchit summary exists for that run-year.") %>%
  tab_footnote("Step 3 uses the INLA model's average marginal effect (accounts for matching covariates), not the raw BACI count directly. It scales the matched-sample AME to the full leakage zone by woody_veg_leakage_zone_ha / (new protected area); per-year bounds are quantiles of delta_area x that positive constant ratio.") %>%
  tab_footnote("Pooled rows: per-year delta_area draws scaled to the leakage zone (Step 3), then within each of 10,000 posterior draws summed across years (Step 3) or combined as sum(leakage)/sum(avoided clearing) (Step 4 aggregate); 2.5/50/97.5% quantiles are then taken over draws. New protected area and avoided clearing on pooled rows are summed over the pooled years. 'Pooled QC (n)' keeps only the n years meeting all 3 matching-QC criteria and is the pooling reported in Figure 5; 'Pooled all (n)' uses every fitted year. The pooled Multiplier is an aggregate ratio, not the pooled median effect divided by a clearing value.")

gtsave(derivation_table, file.path(PLOT_DIR, "figure5_derivation_table.html"))
tryCatch(
  gtsave(derivation_table, file.path(PLOT_DIR, "figure5_derivation_table.png"),
         vwidth = 2100, vheight = 1300),
  error = function(e) message("figure5_derivation_table.png skipped (no headless browser): ",
                              conditionMessage(e))
)

# ── Same table as a publication-ready Excel workbook ─────────────────────────
# Built by hand with {openxlsx}: one worksheet, a merged title/subtitle, merged
# Step / QC spanner headers over row 3, the column headers on row 4, the data
# from row 5; number formats, header fills, a per-cell blue->white->orange
# diverging fill on the Multiplier column, bold-green 95% CI cells where
# significant, frozen panes and footnotes.  Falls back to a plain one-sheet
# dump via {writexl} if {openxlsx} is not installed.
# `Leakage multiplier` is a proportion; blank where avoided clearing is ~0
# (per-year) or n/a (pooled rows).
xlsx_dat <- derivation_df %>%
  dplyr::transmute(
    run = as.character(run), year,
    BC, AC, BI, AI, BACI, BACI_ha,
    ParallelTrends, SMD, PropKeep, QC,
    tte_lb, tte_med, tte_ub, woody_veg_area_ha, five_year_clearing_ha,
    leak_med = if_else(is.finite(leak_med), leak_med, NA_real_),
    sig
  )

xlsx_path <- file.path(PLOT_DIR, "figure5_derivation_table.xlsx")

if (requireNamespace("openxlsx", quietly = TRUE)) {
  ox <- asNamespace("openxlsx")
  n_pl <- sum(xlsx_dat$run == "Leakage to Private Lands")
  n_sf <- sum(xlsx_dat$run == "Leakage to State Forests")
  r0 <- 5L
  r_pl <- r0:(r0 + n_pl - 1L)
  r_sf <- (r0 + n_pl):(r0 + n_pl + n_sf - 1L)
  r_last <- max(r_sf)
  ncol <- 19L

  wb <- ox$createWorkbook()
  ox$addWorksheet(wb, "Leakage multiplier", gridLines = FALSE)

  hdr <- c(
    "Leakage\ndestination", "Year",
    "Before,\ncontrol", "After,\ncontrol", "Before,\ntreated", "After,\ntreated",
    "BACI\n(cells)", "BACI\n(ha)",
    "Parallel\ntrends", "SMD", "Prop.\nkept", "All 3\nmet",
    "2.5%", "Median", "97.5%",
    "New protected\narea (ha woody veg.)", "Avoided\nclearing (ha)", "Leakage\nmultiplier",
    "95% CI"
  )

  ox$writeData(wb, 1, "Derivation of the leakage multiplier", startRow = 1)
  ox$writeData(wb, 1, "From matched-sample clearing counts to the total leakage effect and the multiplier, with matching-QC status and the pooled-over-years estimates", startRow = 2)
  ox$writeData(wb, 1, "Step 1: matched-sample clearing (grid cells)", startCol = 3, startRow = 3)
  ox$writeData(wb, 1, "Step 2: diff-in-diff", startCol = 7, startRow = 3)
  ox$writeData(wb, 1, "Matching QC (model-selected spec)", startCol = 9, startRow = 3)
  ox$writeData(wb, 1, "Step 3: total leakage effect (ha, scaled to full leakage zone)", startCol = 13, startRow = 3)
  ox$writeData(wb, 1, "Step 4: leakage multiplier", startCol = 17, startRow = 3)
  ox$writeData(wb, 1, t(hdr), startRow = 4, colNames = FALSE)
  ox$writeData(wb, 1, xlsx_dat[, -1], startCol = 2, startRow = r0, colNames = FALSE)
  ox$writeData(wb, 1, "Leakage to Private Lands", startCol = 1, startRow = min(r_pl))
  ox$writeData(wb, 1, "Leakage to State Forests", startCol = 1, startRow = min(r_sf))

  ox$mergeCells(wb, 1, cols = 1:ncol, rows = 1)
  ox$mergeCells(wb, 1, cols = 1:ncol, rows = 2)
  ox$mergeCells(wb, 1, cols = 3:6, rows = 3)
  ox$mergeCells(wb, 1, cols = 7:8, rows = 3)
  ox$mergeCells(wb, 1, cols = 9:12, rows = 3)
  ox$mergeCells(wb, 1, cols = 13:16, rows = 3)
  ox$mergeCells(wb, 1, cols = 17:18, rows = 3)
  ox$mergeCells(wb, 1, cols = 1, rows = r_pl)
  ox$mergeCells(wb, 1, cols = 1, rows = r_sf)

  st <- function(...) ox$createStyle(...)
  ox$addStyle(wb, 1, st(fontSize = 14, textDecoration = "bold", halign = "center"), rows = 1, cols = 1:ncol, gridExpand = TRUE)
  ox$addStyle(wb, 1, st(fontSize = 10, textDecoration = "italic", fontColour = "#666666", halign = "center"), rows = 2, cols = 1:ncol, gridExpand = TRUE)
  ox$addStyle(wb, 1, st(fontSize = 9, textDecoration = "bold", halign = "center", valign = "center", fgFill = "#DCE6F1", border = "TopBottomLeftRight", borderColour = "#BFBFBF", wrapText = TRUE), rows = 3, cols = 1:ncol, gridExpand = TRUE)
  ox$addStyle(wb, 1, st(fontSize = 9, textDecoration = "bold", halign = "center", valign = "bottom", fgFill = "#EDEDED", border = "TopBottomLeftRight", borderColour = "#BFBFBF", wrapText = TRUE), rows = 4, cols = 1:ncol, gridExpand = TRUE)
  ox$addStyle(wb, 1, st(fontSize = 10, textDecoration = "bold", valign = "center", fgFill = "#F7F7F7", border = "TopBottomLeftRight", borderColour = "#D9D9D9", wrapText = TRUE), rows = r0:r_last, cols = 1, gridExpand = TRUE)
  ox$addStyle(wb, 1, st(fontSize = 10, halign = "center", border = "TopBottomLeftRight", borderColour = "#E6E6E6"), rows = r0:r_last, cols = 2, gridExpand = TRUE)
  ox$addStyle(wb, 1, st(fontSize = 10, numFmt = "#,##0", halign = "right", border = "TopBottomLeftRight", borderColour = "#E6E6E6"), rows = r0:r_last, cols = c(3:8, 13:17), gridExpand = TRUE)
  ox$addStyle(wb, 1, st(fontSize = 10, halign = "center", border = "TopBottomLeftRight", borderColour = "#E6E6E6"), rows = r0:r_last, cols = 9:12, gridExpand = TRUE)
  ox$addStyle(wb, 1, st(fontSize = 10, halign = "center", border = "TopBottomLeftRight", borderColour = "#E6E6E6"), rows = r0:r_last, cols = 19, gridExpand = TRUE)

  # Diverging fill on the multiplier column (col 18): blue (-200%) -> white (0)
  # -> orange (+200%); grey when off that scale.
  ramp <- grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(201)
  for (k in seq_len(nrow(xlsx_dat))) {
    v <- xlsx_dat$leak_med[k]
    if (is.na(v)) next
    fill <- if (v < -2 || v > 2) "#D9D9D9" else ramp[round((v + 2) / 4 * 200) + 1L]
    fg <- if (abs(v) > 1.1) "#FFFFFF" else "#000000"
    ox$addStyle(wb, 1,
      st(fontSize = 10, numFmt = "0%", halign = "right", fgFill = fill, fontColour = fg,
         border = "TopBottomLeftRight", borderColour = "#E6E6E6"),
      rows = r0 + k - 1L, cols = 18)
  }
  na_rows <- r0 + which(is.na(xlsx_dat$leak_med)) - 1L
  if (length(na_rows)) {
    ox$addStyle(wb, 1, st(fontSize = 10, numFmt = "0%", halign = "right", border = "TopBottomLeftRight", borderColour = "#E6E6E6"), rows = na_rows, cols = 18, gridExpand = TRUE)
  }
  sig_rows <- r0 + which(xlsx_dat$sig == "significant") - 1L
  if (length(sig_rows)) {
    ox$addStyle(wb, 1, st(fontSize = 10, halign = "center", textDecoration = "bold", fontColour = "#1B7837", border = "TopBottomLeftRight", borderColour = "#E6E6E6"), rows = sig_rows, cols = 19, gridExpand = TRUE)
  }
  # Pooled rows in italic (they carry no BC/AC/BI/AI, QC or single denominator).
  pool_rows_idx <- r0 + which(grepl("^Pooled", xlsx_dat$year)) - 1L
  if (length(pool_rows_idx)) {
    ox$addStyle(wb, 1, st(fontSize = 10, textDecoration = "italic", halign = "center", border = "TopBottomLeftRight", borderColour = "#E6E6E6"), rows = pool_rows_idx, cols = 2, gridExpand = TRUE)
  }

  ox$writeData(wb, 1, "BACI = (After, treated - Before, treated) - (After, control - Before, control); a cell is 625 m^2 (25 x 25 m).", startRow = r_last + 2L)
  ox$writeData(wb, 1, "Matching QC: parallel pre-trends (p >= 0.05), SMD (< 0.25), proportion of treated units retained (>= 0.8), on the model-selected matching spec. 'All 3 met' = meets_criteria; blank where no matchit summary exists.", startRow = r_last + 3L)
  ox$writeData(wb, 1, "Step 3 uses the INLA model's matching-covariate-adjusted average marginal effect, not the raw BACI count. Multiplier blank where avoided clearing is negligible; grey fill where it is off the +/-200% scale.", startRow = r_last + 4L)
  ox$writeData(wb, 1, "Pooled rows: per-year delta_area draws scaled to the zone (Step 3), then within each of 10,000 draws summed across years (Step 3) or combined as sum(leakage)/sum(avoided clearing) (Step 4 aggregate); quantiles taken over draws. 'Pooled QC (n)' = the n QC-passing years only = the Figure 5 pooling; 'Pooled all (n)' = every fitted year.", startRow = r_last + 5L)
  ox$addStyle(wb, 1, st(fontSize = 8, textDecoration = "italic", fontColour = "#666666"), rows = (r_last + 2L):(r_last + 5L), cols = 1:ncol, gridExpand = TRUE)

  ox$setColWidths(wb, 1, cols = 1, widths = 20)
  ox$setColWidths(wb, 1, cols = 2, widths = 14)
  ox$setColWidths(wb, 1, cols = 3:8, widths = 10)
  ox$setColWidths(wb, 1, cols = 9:12, widths = 10)
  ox$setColWidths(wb, 1, cols = 13:15, widths = 10)
  ox$setColWidths(wb, 1, cols = 16:17, widths = c(16, 12))
  ox$setColWidths(wb, 1, cols = 18, widths = 12)
  ox$setColWidths(wb, 1, cols = 19, widths = 11)
  ox$setRowHeights(wb, 1, rows = c(3, 4), heights = c(30, 44))
  ox$freezePane(wb, 1, firstActiveRow = r0, firstActiveCol = 3)
  ox$saveWorkbook(wb, xlsx_path, overwrite = TRUE)
} else if (requireNamespace("writexl", quietly = TRUE)) {
  message("openxlsx not installed; writing a plain (unformatted) sheet via writexl.")
  writexl::write_xlsx(setNames(xlsx_dat, c(
    "Leakage destination", "Year",
    "Before, control (cells)", "After, control (cells)",
    "Before, treated (cells)", "After, treated (cells)",
    "BACI (cells)", "BACI (ha)",
    "QC: parallel trends", "QC: SMD", "QC: prop. kept", "QC: all 3 met",
    "Total leakage effect 2.5% (ha)", "Total leakage effect median (ha)",
    "Total leakage effect 97.5% (ha)", "New protected area (ha woody veg.)",
    "Avoided clearing (ha)", "Leakage multiplier", "95% CI"
  )), xlsx_path)
} else {
  message("Neither openxlsx nor writexl installed; ", xlsx_path, " skipped.")
}

# ── Heatmap: leakage multiplier by year x destination ───────────────────────
# Discrete x so only the (unevenly spaced) > 10 ha years get evenly-spaced
# columns; complete() so a year present for one channel but not the other still
# gets a grey, unlabelled cell.
heatmap_plt <- effect_df %>%
  ungroup() %>%
  left_join(dplyr::select(qc_df, run, year, five_year_clearing_ha), by = c("run", "year")) %>%
  filter(five_year_clearing_ha > 10) %>%
  mutate(year = factor(year)) %>%
  tidyr::complete(run, year) %>%
  ggplot(aes(x = year, y = run, fill = leak_med)) +
  geom_tile(color = "grey60", linewidth = 0.4) +
  geom_text(aes(label = if_else(is.na(leak_med), "", label_percent(accuracy = 1)(leak_med)),
                fontface = if_else(sig %in% TRUE, "bold", "plain")), size = 3) +
  scale_x_discrete("Year of state forest revocation", expand = c(0, 0)) +
  scale_fill_gradient2("Leakage\nmultiplier", low = "#0072B2", high = "#D55E00", midpoint = 0,
                       labels = label_percent(), limits = c(-2, 2), oob = squish,
                       na.value = "grey90") +
  labs(y = NULL, title = "Leakage multiplier by year",
       caption = "Bold: 95% CI excludes zero.  Grey: no estimate.  Years with < 10 ha avoided clearing omitted.") +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        plot.caption = element_text(hjust = 0))

ggsave(file.path(PLOT_DIR, "figure5_derivation_heatmap.png"), heatmap_plt,
       width = 1700, height = 750, units = "px", dpi = 300)

message("Derivation table + heatmap saved to plots/")
