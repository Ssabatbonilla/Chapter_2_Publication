# =============================================================================
# 06_figures.R — the data figures, drawn from the tables of this run
# =============================================================================
# WHAT THIS STEP DRAWS (every figure at journal final size, PNG 600 dpi, SVG and PDF, in output/figures/)
# Inputs: the objects earlier steps leave in the same R session (dat, conc, screen, den_mat), the Poor | Fair
# cut-points (H1_BIBI_thresholds_used.csv in the data folder), and the step 11, 15 and 15b tables in output/tables/.
#
#   Main-text figures of this paper (step 18 copies them into output/manuscript/figures/):
#     fig3_hyp12_responsive         Figure 3: every metric of hypotheses 1 and 2 that met the criterion in at least one
#                                   region, as scatters within region, the two weightings of a metric overlaid, with
#                                   each region's Poor | Fair cut drawn
#     fig4_hyp3_responsive          Figure 4: the same for hypothesis 3 (the two weightings, the totals, inequality and
#                                   the two divergences)
#     fig5_effect_sizes_by_region   Figure 5: every screened metric's Spearman rho against the condition gradient with
#                                   its interval, by region and weighting (a, b)
#   Drawn for the companion trait-state manuscript, not this paper (step 18b moves them to trait_state_manuscript/):
#     fig_traitstate_tiers          legacy Figure 6: every trait state associated with condition within its trait group
#                                   under all five seeds, its fourth-corner statistic by region and weighting (step 15)
#     figS4_traitstate_difference   legacy Figure S2: the density-against-biomass difference of each trait state (step 15b)
#   Drawn and not carried into any manuscript:
#     fig4_traitspace_binned        trait space of the sites in the bottom and top third of each region's gradient
#     figS3_structure_scatters, figS4_function_scatters   scatter grids of selected structure and function metrics
#   Working figures read by nothing downstream, drawn only with CH2_WORKING_FIGS=1 (talks and checking):
#     fig3_structure_and_standing_stock, fig5_functional_responses, slide_structure_wide, slide_function_wide
#   Figure S1 of the paper (figS5_chironomid_removal) is drawn by steps 20 and 20b, not here.
#
#   Rules the figures follow: hue encodes region only; condition class is shape, except in Figures 3 and 4 and the
#   trait-state figure, where shape encodes weighting; every metric name, symbol and axis label comes from METRIC_DICT
#   in 00_config.R; no statistic is printed inside a panel; the figures show functional dispersion and not Rao's
#   quadratic entropy, which duplicates it (rho 0.99).
# =============================================================================

say("STEP 06 -- figures")

# CLASS_LEVELS and CLASS_SHAPES now live in 00_config.R so that any step can
# use them; step 19 needs CLASS_LEVELS and cannot source this file.

# theme_pub() in 00_config.R sets journal final size; theme_ch2 adds this step's grid, strip and legend settings.
theme_ch2 <- theme_pub() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey94", colour = NA),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom", legend.box = "horizontal")

plot_dat <- dat %>%
  left_join(conc, by = "SiteID") %>%
  mutate(BIBI_Class = factor(BIBI_Class, levels = CLASS_LEVELS))

# Axis domains derived from the data actually being plotted: 0 to the rounded
# maximum for non-negative metrics, rounded min to max including zero for
# signed ones. Shared by every panel of a metric, across regions and both
# currencies.
Y_DOMAIN <- compute_y_domains(plot_dat)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# THE FIGURE SCHEME (27 Aug 2026)
#   - Panel titles are neutral metric names, not interpretive claims.
#   - Confidence ribbons take the region color (`fill` is mapped).
#   - No figure title inside the graphic; that is the caption's job.
#   - Drawn at journal final size.
#   - Statistics are not drawn on the panels; step 18 generates the coefficients
#     for the captions from the same screen table.
#   - One x-axis title, on the bottom panel only.
#   - Legend stacked vertically: a single horizontal row of region colors and
#     four condition classes is wider than 174 mm and was clipped at both edges.
# -----------------------------------------------------------------------------

# Which regions earned a fitted line. Lines are drawn only where the metric met
# the screening criterion, so a line never implies an association the screen did
# not support. The caption says so.
responsive_regions <- function(metric) {
  screen %>% filter(Metric == metric, Tier == "strong") %>%
    pull(Region) %>% as.character()
}

# ---- The Poor | Fair rating threshold of each region (3 Oct 2026) ---------------------------------------------------
# Drawn on every scatter as a thin gray vertical, as a LABEL of the index's own classes and never as a test:
# Methods 2.5 says classes label sites in figures only. The cut is read from the class file's threshold table
# (H1_BIBI_thresholds_used.csv; Smith et al. 2017, Table 5, family level), never typed.
# Valley and Ridge is scored against two bioregions (29 streams SGV, one SRV); the line is the bioregion that
# scores most of the region's streams, and the caption says so.
# The file comes in two forms: one row per bioregion that scores a stream (Typology, Label, Score, Bioregion,
# n_sites) or one row per region (Typology, Label, Score). Either is read; with the first, the bioregion that scores
# most of a region's streams is kept (Valley and Ridge: SGV, 29 of 30).
THRESH <- read_csv(file.path(H1, "H1_BIBI_thresholds_used.csv"), show_col_types = FALSE) %>% filter(Label == "Poor–Fair")
if (all(c("Bioregion", "n_sites") %in% names(THRESH))) {
  THRESH <- THRESH %>% group_by(Typology) %>% slice_max(n_sites, n = 1, with_ties = FALSE) %>% ungroup()
} else {
  THRESH <- THRESH %>% mutate(Bioregion = NA_character_, n_sites = NA_integer_)
}
THRESH <- THRESH %>% transmute(Region = factor(Typology, levels = REGIONS), cut = Score, Bioregion, n_sites)
stopifnot(nrow(THRESH) == 3, setequal(as.character(THRESH$Region), REGIONS), all(THRESH$cut > 0 & THRESH$cut < 100))
say("  Poor | Fair thresholds drawn: %s", paste(sprintf("%s %.2f", THRESH$Region, THRESH$cut), collapse = "; "))
geom_threshold <- function() {
  geom_vline(data = THRESH, aes(xintercept = cut), colour = "grey45", linetype = "32", linewidth = 0.35)
}

panel_metric <- function(metric, show_x = FALSE, log_y = FALSE) {
  d    <- plot_dat %>% select(Region, BIBI, BIBI_Class, y = all_of(metric))
  resp <- responsive_regions(metric)
  ggplot(d, aes(BIBI, y)) +
    { if (length(resp))
        geom_smooth(data = d %>% filter(Region %in% resp),
                    aes(colour = Region, fill = Region),
                    method = "lm", formula = y ~ x, se = TRUE,
                    linewidth = 0.55, alpha = 0.18, show.legend = FALSE)
      else NULL } +
    geom_point(aes(colour = Region, shape = BIBI_Class), size = 1.35, alpha = 0.9) +
    # NOT scales = "free_y". A free scale overrides the fixed limits set by
    # scale_y_metric(), so each region silently reverted to its own axis and
    # the top Valley & Ridge density point was cut off. A shared 0-to-maximum
    # axis is the point: it shows where each region sits in the full range.
    facet_wrap(~ Region, nrow = 1) +
    scale_colour_manual(values = REGION_COLS, limits = REGIONS,
                        drop = FALSE, guide = "none") +
    scale_fill_manual(values = REGION_COLS, limits = REGIONS,
                      drop = FALSE, guide = "none") +
    scale_shape_manual(values = CLASS_SHAPES, name = "Condition class", drop = FALSE) +
    scale_x_bibi(if (show_x) BIBI_AXIS_LABEL else NULL) +
    scale_y_metric(metric, n = 4) +
    labs(y = mdict_axis(metric), title = wrap_title(mdict_title(metric))) +
    guides(colour = guide_legend(order = 1, nrow = 1, override.aes = list(size = 2.2)),
           shape  = guide_legend(order = 2, nrow = 1, override.aes = list(size = 2))) +
    theme_ch2 +
    theme(plot.title = element_text(size = FIG_BASE_PT - 1.5, face = "plain",
                                    lineheight = 1.05, margin = margin(b = 2)),
          axis.title.x = element_text(size = FIG_BASE_PT - 1.5),
          axis.title.y = element_text(size = FIG_BASE_PT - 1))
}

# ---- Working figure: structure and standing stock (fig3_structure_and_standing_stock) ------
# Four panels; saved only with CH2_WORKING_FIGS=1.
suppressPackageStartupMessages(library(patchwork))
f3 <- (panel_metric("Richness_alpha") /
       panel_metric("HBI_rel_abund") /
       panel_metric("Total_Density") /
       panel_metric("Total_Biomass",  show_x = TRUE)) +
  plot_layout(guides = "collect")
# One legend row only: the region colors are already named by the facet strips.
f3 <- tag_panels(f3) & theme(legend.position = "bottom")
WORKING_FIGS <- nzchar(Sys.getenv("CH2_WORKING_FIGS"))   # 27 Sep 2026: figures nothing reads are opt-in
if (WORKING_FIGS) save_fig(f3, "fig3_structure_and_standing_stock", FIG_W_DOUBLE, 5.6)

# ---- Binned trait space (fig4_traitspace_binned; drawn, not carried into the paper) ------
# The ordination is the one step 11 writes from the current trait matrix and taxon set (11_traitspace_*.csv).
# Step 11 runs before this step.
ss11     <- read_csv(file.path(OUT_TAB, "11_traitspace_site_scores.csv"), show_col_types = FALSE)
pcoa_den <- ss11 %>% filter(Weighting == "Density") %>% select(SiteID, A1, A2)
pcoa_bio <- ss11 %>% filter(Weighting == "Biomass") %>% select(SiteID, A1, A2)
taxa_xy  <- read_csv(file.path(OUT_TAB, "11_traitspace_taxon_coordinates.csv"), show_col_types = FALSE)
axes11   <- read_csv(file.path(OUT_TAB, "11_traitspace_axes.csv"), show_col_types = FALSE)
stopifnot(setequal(taxa_xy$TAXANAME, colnames(den_mat)))

# The site scores carry no B-IBI; the canonical B-IBI is joined from `dat`.
prep_pcoa <- function(x, weighting) {
  x %>% mutate(SiteID = std_site(SiteID)) %>%
    select(SiteID, A1, A2) %>%
    left_join(dat %>% select(SiteID, Region, BIBI), by = "SiteID") %>%
    mutate(Weighting = weighting)
}
pc <- bind_rows(prep_pcoa(pcoa_den, "Density"), prep_pcoa(pcoa_bio, "Biomass")) %>%
  mutate(Region = factor(as.character(Region), levels = REGIONS)) %>%
  droplevels() %>%
  group_by(Region) %>%
  mutate(Condition = case_when(
    BIBI <= quantile(BIBI, 1/3) ~ "Low condition (bottom third)",
    BIBI >= quantile(BIBI, 2/3) ~ "High condition (top third)",
    TRUE ~ NA_character_)) %>%
  ungroup() %>%
  filter(!is.na(Condition)) %>%
  mutate(Condition = factor(Condition,
           levels = c("Low condition (bottom third)", "High condition (top third)")),
         Weighting = factor(Weighting, levels = c("Density", "Biomass")))

hulls <- pc %>% group_by(Region, Weighting, Condition) %>%
  slice(chull(A1, A2)) %>% ungroup()
cents <- pc %>% group_by(Region, Weighting, Condition) %>%
  summarise(A1 = mean(A1), A2 = mean(A2), .groups = "drop")
arrows_df <- cents %>%
  pivot_wider(names_from = Condition, values_from = c(A1, A2)) %>%
  rename(x = 3, y = 5, xend = 4, yend = 6)

# fdiv_lab is not used in any figure this step saves.
fdiv_lab <- screen %>% filter(Metric == "FDiv_den") %>%
  transmute(Region = as.character(Region), Weighting = factor("Density", levels=c("Density","Biomass")),
            lab = sprintf("FDiv (density) rho = %+.2f, FDR = %s", rho,
                          ifelse(p_FDR < 0.001, "<0.001", sprintf("%.3f", p_FDR))))

# POTENTIAL VERSUS REALIZED POOL.
#   Each region holds a different share of the metacommunity's taxa (repository table
#   regional_pools_shared_taxa.csv), so the backdrop shows which part of trait space a region can occupy.
#   Open circles are taxa in the metacommunity that were never collected in
#   that region (potential); filled circles are taxa present there (realized).
#   Shape is used rather than color because color already encodes condition
#   class, and one aesthetic must mean one thing across the figure set.
taxa_pool <- purrr::map_dfr(REGIONS, function(rg) {
  s    <- dat$SiteID[dat$Region == rg]
  pres <- colnames(den_mat)[colSums(den_mat[s, , drop = FALSE]) > 0]
  taxa_xy %>%
    mutate(Region = factor(rg, levels = REGIONS),
           Pool = ifelse(TAXANAME %in% pres,
                         "Present in this region (realized)",
                         "Metacommunity only (not collected here)"))
}) %>%
  mutate(Region = droplevels(Region)) %>%
  tidyr::crossing(Weighting = factor(c("Density", "Biomass"),
                                     levels = c("Density", "Biomass"))) %>%
  mutate(Pool = factor(Pool, levels = c("Metacommunity only (not collected here)",
                                        "Present in this region (realized)")))

f4 <- ggplot(pc, aes(A1, A2)) +
  # Three aesthetics differ at once -- shape, size and alpha -- because a
  # single one was not enough: open and filled circles at the same gray read as
  # one cloud. Potential is a faint small ring, realized a solid dark dot.
  geom_point(data = taxa_pool,
             aes(A1, A2, shape = Pool, size = Pool, alpha = Pool),
             colour = "grey20", inherit.aes = FALSE) +
  scale_shape_manual(values = c(1, 16), name = NULL) +
  scale_size_manual(values = c(0.5, 1.15), name = NULL) +
  scale_alpha_manual(values = c(0.30, 0.95), name = NULL) +
  geom_polygon(data = hulls, aes(fill = Condition, colour = Condition, group = Condition),
               alpha = 0.20, linewidth = 0.5) +
  geom_point(aes(colour = Condition), size = 1.6) +
  geom_segment(data = arrows_df, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.16, "cm")), linewidth = 0.6,
               colour = "black", inherit.aes = FALSE) +
  facet_grid(Region ~ Weighting) +
  scale_colour_manual(values = PAL_CONDITION) +
  scale_fill_manual(values = PAL_CONDITION) +
  # No baked-in caption and no in-panel coefficients: both move to the generated
  # caption in step 18, where they can be edited without redrawing the figure.
  scale_x_continuous(labels = NEG_LAB) + scale_y_continuous(labels = NEG_LAB) +
  labs(x = sprintf("PCoA axis 1 (%.1f%%)", axes11$pct_variance[1]), y = sprintf("PCoA axis 2 (%.1f%%)", axes11$pct_variance[2]),
       colour = NULL, fill = NULL) +
  guides(shape = guide_legend(order = 1,
           override.aes = list(size = c(1.5, 2.1), alpha = c(0.45, 0.95))),
         size = "none", alpha = "none",
         colour = guide_legend(order = 2), fill = guide_legend(order = 2)) +
  theme_ch2 + theme(legend.position = "bottom", legend.box = "vertical",
                    legend.spacing.y = grid::unit(0.02, "cm"))
# Single composite panel -- no (a) tag; it would label nothing.
save_fig(f4, "fig4_traitspace_binned", FIG_W_DOUBLE, 5.4)

# ---- Working figure: functional responses (fig5_functional_responses; CH2_WORKING_FIGS=1 only) ----
# Same scheme as panel_metric: neutral title, no in-panel statistics, region
# color on both the line and its ribbon. The x-axis title appears once, on the
# bottom row only.
resp_panel <- function(metric, show_x = FALSE, log_y = FALSE) {
  d    <- plot_dat %>% select(Region, BIBI, BIBI_Class, y = all_of(metric))
  resp <- responsive_regions(metric)
  ggplot(d, aes(BIBI, y)) +
    { if (length(resp))
        geom_smooth(data = d %>% filter(Region %in% resp),
                    aes(colour = Region, fill = Region),
                    method = "lm", formula = y ~ x, se = TRUE,
                    linewidth = 0.55, alpha = 0.18, show.legend = FALSE)
      else NULL } +
    geom_point(aes(colour = Region, shape = BIBI_Class), size = 1.3, alpha = 0.9) +
    # These panels have no region facets, so color is the ONLY region encoding here
    # and needs a key. panel_metric() is faceted by region and suppresses it.
    scale_colour_manual(values = REGION_COLS, limits = REGIONS, drop = FALSE,
                        name = "Region") +
    scale_fill_manual(values = REGION_COLS, limits = REGIONS, drop = FALSE,
                      guide = "none") +
    scale_shape_manual(values = CLASS_SHAPES, name = "Condition class", drop = FALSE) +
    scale_x_bibi(if (show_x) BIBI_AXIS_LABEL else NULL) +
    scale_y_metric(metric, n = 4) +
    labs(y = mdict_axis(metric), title = wrap_title(mdict_title(metric))) +
    guides(shape = guide_legend(nrow = 1, override.aes = list(size = 2))) +
    theme_ch2 +
    theme(plot.title = element_text(size = FIG_BASE_PT - 1.5, face = "plain",
                                    lineheight = 1.05, margin = margin(b = 2)),
          axis.title.x = element_text(size = FIG_BASE_PT - 1.5),
          axis.title.y = element_text(size = FIG_BASE_PT - 1))
}

# Titles name the METRIC and the property it belongs to. They do not state the
# result -- that is the caption's job.
# PANEL (b) IS FDis, NOT RaoQ. The two correlate at rho = 0.99 in these data
# (density and biomass alike), so they are one signal, not two. Methods 2.4
# names functional dispersion the primary measure of assemblage spread and
# Rao's Q the concordance check; the figure must show the primary one.
f5 <- (resp_panel("FDiv_den")  + resp_panel("FDis_den")) /
      (resp_panel("FR_RD_den") + resp_panel("FR_RD_bio")) /
      (resp_panel("FR_DB_RD", show_x = TRUE) +
       resp_panel("Gini_DB_DOM", show_x = TRUE)) +
  plot_layout(guides = "collect")
f5 <- tag_panels(f5) & theme(legend.position = "bottom", legend.box = "vertical",
                             legend.box.just = "center")
if (WORKING_FIGS) save_fig(f5, "fig5_functional_responses", FIG_W_DOUBLE, 6.8)

# ---- resp_grid(): one row per metric, one column per region ----------------------------
# Used for the scatter grids figS3_structure_scatters and figS4_function_scatters (end of this step).
# REGION FACETS, NOT A REGION OVERLAY (2 September 2026).
#   Overlaying three regions in one panel puts them on one pair of axes, which
#   reads as a cross-region comparison. This paper never runs one: every test is
#   within region, because the index is calibrated bioregion by bioregion.
#   Faceting by region says so -- three panels, three analyses -- and each row
#   has its own y scale (scales = "free_y"), so the three regions of one metric
#   share a scale while never sharing an axis.
#
#   Region color is kept because the palette is document-wide, but its legend is
#   suppressed: the strip already names the region.
resp_grid <- function(metrics, height_in) {
  miss <- setdiff(metrics, names(ROW_LABELS))
  if (length(miss))
    stop("resp_grid: no ROW_LABELS entry for ", paste(miss, collapse = ", "),
         ". A missing label becomes NA and silently merges facets. Add it.")
  long <- purrr::map_dfr(metrics, function(m) {
    plot_dat %>%
      transmute(Region, BIBI, BIBI_Class, metric = m, y = .data[[m]])
  }) %>%
    mutate(metric = factor(metric, levels = metrics,
                           labels = ROW_LABELS[metrics]))

  fit <- purrr::map_dfr(metrics, function(m) {
    rg <- responsive_regions(m)
    if (!length(rg)) return(NULL)
    plot_dat %>% filter(Region %in% rg) %>%
      transmute(Region, BIBI, metric = m, y = .data[[m]])
  })
  if (nrow(fit))
    fit <- fit %>% mutate(metric = factor(metric, levels = metrics,
                          labels = ROW_LABELS[metrics]))

  ggplot(long, aes(BIBI, y)) +
    geom_threshold() +
    { if (nrow(fit))
        geom_smooth(data = fit, aes(colour = Region, fill = Region),
                    method = "lm", formula = y ~ x, se = TRUE,
                    linewidth = 0.5, alpha = 0.18, show.legend = FALSE)
      else NULL } +
    geom_point(aes(colour = Region, shape = BIBI_Class), size = 0.95, alpha = 0.9) +
    facet_grid(metric ~ Region, scales = "free_y", switch = "y") +
    scale_colour_manual(values = REGION_COLS, limits = REGIONS, drop = FALSE,
                        guide = "none") +
    scale_fill_manual(values = REGION_COLS, limits = REGIONS, drop = FALSE,
                      guide = "none") +
    scale_shape_manual(values = CLASS_SHAPES, name = "Condition class", drop = FALSE) +
    scale_x_bibi(BIBI_AXIS_LABEL) +
    scale_y_continuous(n.breaks = 3, labels = NEG_LAB) +
    labs(y = NULL) +
    guides(shape = guide_legend(nrow = 1, override.aes = list(size = 2))) +
    theme_ch2 +
    theme(legend.position = "bottom",
          strip.placement  = "outside",
          strip.text.y.left = element_text(angle = 90, size = FIG_BASE_PT - 2.5,
                                           lineheight = 0.95),
          strip.text.x     = element_text(size = FIG_BASE_PT - 1, face = "plain"),
          panel.spacing.x  = unit(3, "pt"),
          panel.spacing.y  = unit(3.5, "pt"),
          axis.text        = element_text(size = FIG_BASE_PT - 2.5))
}

# Row strips carry the metric and its currency on two lines. Units are given in
# the caption; a strip is 20 mm and cannot hold "mg DM m-2".
ROW_LABELS <- c(
  Richness_alpha = "Taxonomic\nrichness",
  FRic_den       = "Functional\nrichness",
  Shannon        = "Shannon\ndiversity (den.)",
  HBI_rel_abund  = "Tolerance\n(rel. abundance)",
  HBI_bio        = "Tolerance\n(biomass)",
  Total_Density  = "Total density\n(density)",
  Total_Biomass  = "Total biomass\n(biomass)",
  FDiv_den       = "Functional\ndivergence (den.)",
  FDis_den       = "Functional\ndispersion (den.)",
  FR_RD_den      = "Overlap\nredundancy (den.)",
  Gini_DB_DOM    = "Concentration\ndivergence",
  FR_DB_RD       = "Redundancy\ndivergence",
  FR_RD_bio      = "Overlap\nredundancy (bio.)"
)

# ---- Figure 5: the whole screen as one compact effect-size figure (fig5_effect_sizes_by_region) ----
#   Every metric-by-weighting test in every region on one axis -- the correlation
#   itself -- so the result is one shape rather than a wall of scatter panels.
#
#   Regions stay faceted, so no cross-region test is implied. The criterion is
#   drawn as a pair of reference lines, and a filled point means the test met it.
#   Confidence intervals are Fisher-z on the Spearman coefficient. They show the
#   precision of each estimate and NOTHING MORE -- do not write that they bound
#   the null. At n = 30 they are wide; many Valley & Ridge intervals still
#   contain the minimum detectable effect. The detection bound is the power calculation
#   (step 14), and the regional null rests on three independent methods
#   agreeing, not on any single interval.
#
#   (a) structure and totals, (b) functional diversity, redundancy, inequality and the divergences.
# HBI_rel_abund is labeled explicitly here; every other metric takes its axis label from the metric
# dictionary (mdict_axis), falling back to the metric key.
axis_lab <- function(m) {
  if (identical(m, "HBI_rel_abund")) return("HBI (rel. abundance)")
  tryCatch(mdict_axis(m), error = function(e) m)
}

fisher_ci <- function(rho, n, conf = 0.95) {
  z  <- atanh(pmax(pmin(rho, 0.999), -0.999))
  se <- 1.06 / sqrt(n - 3)                       # the Fieller, Hartley and Pearson (1957) approximation for Spearman's rho;
                                                 # Bonett and Wright (2000) give sqrt((1 + rho^2/2)/(n - 3)), which equals
                                                 # this only at |rho| = 0.5
  k  <- qnorm(1 - (1 - conf) / 2)
  list(lo = tanh(z - k * se), hi = tanh(z + k * se))
}

BLOCK_A <- c("Richness_alpha","Shannon","HBI_rel_abund","HBI_bio",
             "Total_Density","Total_Biomass")
BLOCK_B <- c("FRic_den","FEve_den","FEve_bio","FDiv_den","FDiv_bio",
             "FDis_den","FDis_bio","FR_RD_den","FR_RD_bio","FR_SF_den_log10",
             "Gini_den","Gini_bio","N50_den","N50_bio","Gini_DB_DOM","FR_DB_RD")

eff <- screen %>%
  filter(Metric %in% c(BLOCK_A, BLOCK_B)) %>%
  mutate(ci = list(NULL)) %>% rowwise() %>%
  mutate(lo = fisher_ci(rho, n)$lo, hi = fisher_ci(rho, n)$hi) %>% ungroup() %>%
  mutate(
    Block = ifelse(Metric %in% BLOCK_A,
                   "(a) Structure\nand totals",
                   "(b) Functional diversity, redundancy,\ninequality and the divergences"),
    Label = vapply(Metric, axis_lab, character(1)),
    Label = factor(Label, levels = rev(vapply(c(BLOCK_A, BLOCK_B), axis_lab, character(1)))),
    Met   = Tier == "strong")

f3 <- ggplot(eff, aes(rho, Label, colour = Region)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_vline(xintercept = c(-RHO_MIN, RHO_MIN), colour = "grey72",
             linetype = "22", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.42, alpha = 0.75) +
  geom_point(aes(shape = Met, fill = Region), size = 1.5, stroke = 0.45) +
  facet_grid(Block ~ Region, scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = REGION_COLS, limits = REGIONS, guide = "none") +
  scale_fill_manual(values   = REGION_COLS, limits = REGIONS, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), guide = "none") +
  scale_x_continuous("Spearman correlation with biological condition",
                     limits = c(-1, 1), breaks = seq(-1, 1, 0.5), labels = NEG_LAB) +
  labs(y = NULL) +
  theme_ch2 +
  theme(strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 90, size = FIG_BASE_PT - 2.5),
        strip.text.x      = element_text(size = FIG_BASE_PT - 1),
        panel.spacing.x   = unit(3, "pt"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25),
        axis.text.y       = element_text(size = FIG_BASE_PT - 2.5),
        axis.text.x       = element_text(size = FIG_BASE_PT - 2))

# ---- Trait-state difference panel (figS4_traitstate_difference; step 15b) ----------------------------------------------
# Legacy Figure S2; it belongs to the companion trait-state manuscript, not this paper. One row per counted contrast
# (one state per two-state group, 38 per region), the difference weighted by biomass minus weighted by density, the 95%
# percentile interval as its range over the five seeds, an open diamond where the state is not shown to differ; a filled
# diamond is reserved for a difference that passes under every seed. Trait-state names as Twardochleb et al. publish them.
wdiff <- read_csv(file.path(OUT_TAB, "15b_weighting_difference.csv"), show_col_types = FALSE) %>%
  filter(keep) %>% mutate(code = sub("^BIBI / ", "", term))
stopifnot(all(wdiff$code %in% names(TRAIT_STATE_LAB)), all(wdiff$block %in% TRAIT_GROUP_ORDER))
wdiff <- wdiff %>%
  mutate(Group = factor(TRAIT_GROUP_SHORT[block], levels = TRAIT_GROUP_SHORT[TRAIT_GROUP_ORDER]),
         State = TRAIT_STATE_LAB[code], Region = factor(Region, levels = REGIONS),
         Differs = verdict == "differs under every seed") %>%
  arrange(Group, match(code, names(TRAIT_STATE_LAB))) %>%
  mutate(State = factor(State, levels = rev(unique(State))))
stopifnot(nrow(wdiff) == 3 * 38)
say("  trait-state difference (figS4_traitstate_difference): %d contrasts, %d differ under every seed", nrow(wdiff), sum(wdiff$Differs))
f3c <- ggplot(wdiff, aes(diff, State, colour = Region)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo95_min, xmax = hi95_max), height = 0, linewidth = 1.6, colour = "grey72", alpha = 0.9) +
  geom_point(aes(shape = Differs, fill = Region), size = 1.6, stroke = 0.45) +
  facet_grid(Group ~ Region, scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = REGION_COLS, limits = REGIONS, guide = "none") +
  scale_fill_manual(values = REGION_COLS, limits = REGIONS, guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 5), guide = "none") +
  scale_x_continuous("Fourth-corner statistic weighted by biomass minus weighted by density",
                     limits = c(-1, 1), breaks = seq(-1, 1, 0.5), labels = NEG_LAB) +
  labs(y = NULL) + theme_ch2 +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = FIG_BASE_PT - 2.5, face = "plain", colour = "grey35"),
        strip.background.y = element_blank(), strip.text.x = element_blank(),
        panel.spacing.x = unit(3, "pt"), panel.spacing.y = unit(1.5, "pt"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25),
        axis.text.y = element_text(size = FIG_BASE_PT - 2.5), axis.text.x = element_text(size = FIG_BASE_PT - 2))
# Panels (a) and (b) are Figure 5; the trait-state difference panel is saved on its own. The unlink() calls remove
# the files of retired stems (fig3_ and fig4_effect_sizes_by_region) so no stale copy survives.
save_fig(f3, "fig5_effect_sizes_by_region", FIG_W_DOUBLE, 5.0)
unlink(file.path(OUT_FIG, paste0("fig4_effect_sizes_by_region", c(".png", ".pdf", ".svg"))))
save_fig(f3c + theme(strip.text.x = element_text(size = FIG_BASE_PT - 1, face = "bold")),
         "figS4_traitstate_difference", FIG_W_DOUBLE, 6.6)
unlink(file.path(OUT_FIG, paste0("fig3_effect_sizes_by_region", c(".png", ".pdf", ".svg"))))   # no stale three-panel copy survives

# ---- Figures 3 and 4: the responsive metrics, by hypothesis (3 Oct 2026) ---------------------------------------------
# The data figures show every metric that responded, organized by the hypothesis that tests it. THE RULE, asserted
# below: a metric (or a density-biomass pair of the same metric) is drawn if it met the screening criterion in at least
# one region under at least one weighting; a pair is drawn whole so the weighting that did not respond is seen beside
# the one that did; and every metric that met the criterion anywhere must land in some row (Rao's quadratic entropy
# excepted: it duplicates functional dispersion, rho 0.99, and the caption says so). Figure 3 is hypotheses 1 and 2
# (structure; traits against structure), Figure 4 is hypothesis 3 (the two weightings, the totals, inequality and the
# two divergences). Encoding: hue is region and nothing else; within a row the two weightings are filled/solid
# (density) against open/dashed (biomass); a fitted line is drawn only where that metric and weighting met the
# criterion; no statistic inside a panel; each region's Poor | Fair cut is a gray vertical label (THRESH).
# Condition-class shapes are not drawn in these two figures because shape encodes weighting; the caption says so.
# Rows have free y scales; the two totals are separate rows because their units differ.
HYP_ROWS <- list(
  fig3 = list(
    list(label = "Taxonomic\nrichness",      den = "Richness_alpha",  bio = NA),
    list(label = "Shannon\ndiversity",       den = "Shannon",         bio = NA),
    list(label = "Pollution\ntolerance (HBI)", den = "HBI_rel_abund", bio = "HBI_bio"),
    list(label = "Functional\nrichness",     den = "FRic_den",        bio = NA),
    list(label = "Space-filling\nredundancy (log10)", den = "FR_SF_den_log10", bio = NA),
    list(label = "Functional\ndispersion",   den = "FDis_den",        bio = "FDis_bio")),
    # Functional divergence (FDiv_den, FDiv_bio) was a Figure 3 row until 6 October 2026. Under the density method of
    # 5 October (James et al. 2022: mean over processed replicates, no occurrence weighting) it meets the criterion in
    # no region under either weighting (Valley and Ridge FDiv_bio FDR 0.042 -> 0.218), so the rule below removes it.
  fig4 = list(
    list(label = "Overlap\nredundancy",      den = "FR_RD_den",       bio = "FR_RD_bio"),
    list(label = "Taxa making up\nhalf the total (N50)", den = "N50_den", bio = "N50_bio"),
    list(label = "Gini\ncoefficient",        den = "Gini_den",        bio = "Gini_bio"),
    list(label = "Total density\n(ind. m⁻²)",  den = "Total_Density",  bio = NA, pair = "Total_Biomass"),   # the two totals are one pair (hypothesis 3) in two rows: their units differ
    list(label = "Total biomass\n(mg DM m⁻²)", den = NA,  bio = "Total_Biomass", pair = "Total_Density"),
    list(label = "Redundancy\ndivergence",   den = "FR_DB_RD",        bio = NA),
    list(label = "Concentration\ndivergence", den = "Gini_DB_DOM",   bio = NA))   # 6 Oct 2026: renamed from "Inequality divergence"
)
strong_any <- screen %>% filter(Tier == "strong") %>% pull(Metric) %>% unique()
row_metrics <- function(r) na.omit(c(r$den, r$bio))
for (fg in names(HYP_ROWS)) for (r in HYP_ROWS[[fg]]) {
  stopifnot(all(row_metrics(r) %in% names(plot_dat)), all(row_metrics(r) %in% screen$Metric))
  if (!any(c(row_metrics(r), r$pair) %in% strong_any))   # a row stands if it or its declared pair member responded
    stop(sprintf("%s row '%s' met the criterion in no region under either weighting; the rule admits responsive metrics only", fg, gsub("\n", " ", r$label)))
}
drawn <- unlist(lapply(unlist(HYP_ROWS, recursive = FALSE), row_metrics))
left_out <- setdiff(strong_any, c(drawn, "RaoQ_den", "RaoQ_bio"))
if (length(left_out)) stop("responsive metric(s) in no row of Figure 3 or 4: ", paste(left_out, collapse = ", "))
say("  Figures 3 and 4 rows: %d + %d; every responsive metric drawn (Rao's Q omitted as the duplicate of FDis)",
    length(HYP_ROWS$fig3), length(HYP_ROWS$fig4))

WT_LEVELS <- c("weighted by density", "weighted by biomass")
# 6 Oct 2026: tick labels formatted per panel to three significant digits. NEG_LAB (label_number, accuracy NULL) chose one
# accuracy for a whole free-scaled figure, so Figure 4's Gini row printed "1" at 0.8 and its N50 row lost two of three labels.
LAB_ROW <- function(x) {
  out <- formatC(signif(x, 3), format = "fg", digits = 3, big.mark = " ", drop0trailing = TRUE)
  out <- trimws(sub("\\.$", "", out)); out <- gsub("-", "\u2212", out); out[is.na(x)] <- NA; out
}
# 6 Oct 2026: every row carries the same format, three or four values spanning its range (low, middle, high).
BRK_ROW <- function(lims) {
  spans0 <- lims[1] < 0 && lims[2] > 0
  inside <- function(b) b[b >= lims[1] - 1e-9 & b <= lims[2] + 1e-9]
  for (n in 3:12) {
    b <- inside(pretty(lims, n = n))
    while (length(b) > 4) {                               # thin to every other value, keeping 0 where the row spans it
      z <- if (spans0) which.min(abs(b)) else 1
      b <- b[seq_along(b) %% 2 == z %% 2]
    }
    if (length(b) >= 3 && (!spans0 || any(abs(b) < 1e-12))) return(b)
  }
  r <- function(v) { d <- 10^(floor(log10(abs(v)))); trunc(v / d) * d }   # one significant digit, rounded toward zero (inside)
  if (spans0) return(c(r(lims[1] * 0.9), 0, r(lims[2] * 0.9)))
  inside(signif(seq(lims[1], lims[2], length.out = 3), 2))
}
REGION_STRIP <- c("Valley & Ridge" = "Valley and Ridge", "Piedmont" = "Piedmont", "Coastal Plain" = "Coastal Plain")   # 6 Oct 2026: region names spelled out in the strips
hyp_grid <- function(rows) {
  labs_ <- vapply(rows, `[[`, character(1), "label")
  long <- purrr::map_dfr(rows, function(r) {
    purrr::map_dfr(c("den", "bio"), function(w) {
      m <- r[[w]]; if (is.na(m)) return(NULL)
      plot_dat %>% transmute(Region, BIBI, row = r$label, metric = m, y = .data[[m]],
                             Weighting = if (w == "den") WT_LEVELS[1] else WT_LEVELS[2])
    })
  }) %>% mutate(row = factor(row, levels = labs_), Weighting = factor(Weighting, levels = WT_LEVELS))
  fit <- long %>% group_by(metric) %>% filter(Region %in% responsive_regions(first(metric))) %>% ungroup()
  # relative-abundance HBI is the index's own weighting, drawn as the density member of its pair
  ggplot(long, aes(BIBI, y)) +
    geom_threshold() +
    # The screen is a rank correlation, so the drawn trend must not claim a straight line: a GAM smooth
    # (mgcv, s(BIBI, k = 3), REML), fitted only where that metric and weighting met the criterion.
    { if (nrow(fit)) geom_smooth(data = fit, aes(colour = Region, fill = Region, linetype = Weighting, group = interaction(metric, Region)),
                                 method = "gam", formula = y ~ s(x, k = 3), method.args = list(method = "REML"),
                                 se = TRUE, linewidth = 0.5, alpha = 0.14, show.legend = FALSE) else NULL } +
    geom_point(aes(colour = Region, fill = Region, shape = Weighting), size = 0.95, stroke = 0.4, alpha = 0.9) +
    facet_grid(row ~ Region, scales = "free_y", switch = "y", labeller = labeller(Region = REGION_STRIP)) +
    scale_colour_manual(values = REGION_COLS, limits = REGIONS, drop = FALSE, guide = "none") +
    scale_fill_manual(values = REGION_COLS, limits = REGIONS, drop = FALSE, guide = "none") +
    scale_shape_manual(values = setNames(c(21, 1), WT_LEVELS), name = NULL, drop = FALSE,
                       guide = guide_legend(override.aes = list(fill = c("grey30", NA), colour = "grey30", size = 2))) +
    scale_linetype_manual(values = setNames(c("solid", "22"), WT_LEVELS), guide = "none") +
    scale_x_bibi(BIBI_AXIS_LABEL) +
    # 6 Oct 2026: BRK_ROW and LAB_ROW give each free-scaled row three or four readable ticks; n.breaks = 3 had
    # left Figure 4's N50 and Gini rows with one tick each.
    scale_y_continuous(breaks = BRK_ROW, labels = LAB_ROW) +
    labs(y = NULL) + theme_ch2 +
    theme(legend.position = "bottom", legend.text = element_text(size = FIG_BASE_PT - 2), legend.key.height = unit(8, "pt"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90, size = FIG_BASE_PT - 2.5, lineheight = 0.95),
          strip.text.x = element_text(size = FIG_BASE_PT - 1, face = "plain"),
          panel.spacing.x = unit(3, "pt"), panel.spacing.y = unit(3, "pt"),
          axis.text = element_text(size = FIG_BASE_PT - 2.5))
}
save_fig(hyp_grid(HYP_ROWS$fig3), "fig3_hyp12_responsive", FIG_W_DOUBLE, 8.6)
save_fig(hyp_grid(HYP_ROWS$fig4), "fig4_hyp3_responsive",  FIG_W_DOUBLE, 8.6)
unlink(file.path(OUT_FIG, paste0("fig3_two_pairs", c(".png", ".pdf", ".svg"))))

# ---- Trait-state figure (fig_traitstate_tiers; legacy Figure 6, redrawn 6 Oct 2026) ------------------------------------
# Not in this paper: it belongs to the companion trait-state manuscript (step 18b moves it there). One concept: which
# trait states were associated with the condition gradient, in which direction, under which weighting, in each region.
# Encoded once each: row = trait state (only states passing within their trait group under all five seeds; a two-state
# group shows one state, its mirror is the opposite statistic), column = region (hue repeats region, as in every
# figure), horizontal position = the fourth-corner statistic (sign = direction), filled circle = weighted by density,
# open triangle = weighted by biomass. Left to the tables and the caption: borderline states (some seeds only), the
# region-wide correction, the per-group test, single-taxon dependence. Asserted: every state of the per-block
# fourth-corner tables (legacy Table 4) that passes under every seed is drawn or mirrored.
PASS_FC <- "passes under every seed"
rd6 <- function(f) read_csv(file.path(OUT_TAB, f), show_col_types = FALSE)
pb6 <- bind_rows(mutate(rd6("15_fourthcorner_per_block.csv"), Weighting = "weighted by density"),
                 mutate(rd6("15_fourthcorner_per_block_biomass.csv"), Weighting = "weighted by biomass")) %>%
  mutate(code = sub("^BIBI / ", "", term))
stopifnot(all(pb6$code %in% names(TRAIT_STATE_LAB)), all(pb6$block %in% TRAIT_GROUP_ORDER))
two_state <- pb6 %>% distinct(block, code) %>% count(block) %>% filter(n == 2) %>% pull(block)
pass6 <- pb6 %>% filter(verdict == PASS_FC)
# the state a two-state group is shown by: the one whose passing statistics sum above zero (e.g. weak flight, low dispersal)
canon <- pass6 %>% filter(block %in% two_state) %>% group_by(block, code) %>% summarise(s = sum(stat), .groups = "drop") %>%
  group_by(block) %>% slice_max(s, n = 1, with_ties = FALSE) %>% pull(code)
show6 <- pass6 %>% filter(!(block %in% two_state) | code %in% canon)
mirrored <- pass6 %>% filter(block %in% two_state, !(code %in% canon))
stopifnot(all(paste(mirrored$Region, mirrored$Weighting, mirrored$block) %in% paste(show6$Region, show6$Weighting, show6$block)))
show6 <- show6 %>%
  mutate(row = paste0(TRAIT_GROUP_SHORT[block], ": ", TRAIT_STATE_LAB[code]),
         ord = match(block, TRAIT_GROUP_ORDER) * 100 + match(code, names(TRAIT_STATE_LAB)) / 100,
         Region = factor(Region, levels = REGIONS), Weighting = factor(Weighting, levels = WT_LEVELS))
rows6 <- show6 %>% distinct(row, ord) %>% arrange(desc(ord)) %>% pull(row)
show6 <- show6 %>% mutate(y = match(row, rows6)) %>%
  group_by(Region, row) %>% mutate(dy = if (n() == 2) ifelse(Weighting == WT_LEVELS[1], 0.16, -0.16) else 0) %>% ungroup()
say("  trait-state figure (legacy Figure 6): %d trait states drawn in %d rows (%d mirrored states of two-state groups shown once); %d passing rows of the per-block fourth-corner tables",
    nrow(show6), length(rows6), nrow(mirrored), nrow(pass6))
xlim6 <- max(abs(show6$stat)) * 1.08
f6 <- ggplot(show6, aes(stat, y + dy)) +
  geom_hline(yintercept = seq_along(rows6), colour = "grey92", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_point(aes(colour = Region, fill = Region, shape = Weighting), size = 1.9, stroke = 0.5) +
  facet_grid(. ~ Region, labeller = labeller(Region = REGION_STRIP)) +
  scale_shape_manual(values = setNames(c(21, 2), WT_LEVELS), name = NULL, drop = FALSE,
                     guide = guide_legend(override.aes = list(fill = c("grey30", NA), colour = "grey30", size = 2))) +
  scale_colour_manual(values = REGION_COLS, limits = REGIONS, guide = "none") +
  scale_fill_manual(values = REGION_COLS, limits = REGIONS, guide = "none") +
  scale_x_continuous(limits = c(-xlim6, xlim6), breaks = BRK_ROW,
                     labels = LAB_ROW) +
  scale_y_continuous(breaks = seq_along(rows6), labels = rows6, expand = expansion(add = 0.6)) +
  labs(x = "Association with biological condition (fourth-corner statistic)\nbelow 0: more in worse-condition streams; above 0: more in better-condition streams",
       y = NULL) +
  theme_ch2 +
  theme(legend.position = "bottom", legend.text = element_text(size = FIG_BASE_PT - 2), legend.key.height = unit(8, "pt"),
        panel.grid = element_blank(), strip.text.x = element_text(size = FIG_BASE_PT - 1, face = "plain"),
        axis.text.y = element_text(size = FIG_BASE_PT - 2), axis.text.x = element_text(size = FIG_BASE_PT - 2.5),
        axis.title.x = element_text(size = FIG_BASE_PT - 2, lineheight = 1.05), panel.spacing.x = unit(6, "pt"),
        legend.margin = margin(t = 4))
save_fig(f6, "fig_traitstate_tiers", FIG_W_DOUBLE, 0.95 + 0.2 * length(rows6))

# ---- Scatter grids and slide variants (not carried into the paper) ------------------------------------------------
# figS3_structure_scatters and figS4_function_scatters (saved at the end of this step) show the points behind
# selected coefficients of Figure 5, with the same region facets and the same rule that a fitted line appears only
# where the screen supported one. Shannon diversity is in the structure grid because prediction 1a names it.
# ---- Slide variants: the same grids, transposed ------------------------------
# A 6-row x 3-column grid is portrait; a slide is landscape, and shrinking the
# portrait version to fit makes the panels unreadable. These put regions on rows
# and metrics on columns, which is the same information at a slide's aspect.
# Manuscript figures are unaffected -- these are extra files.
resp_grid_wide <- function(metrics, h_in) {
  miss <- setdiff(metrics, names(ROW_LABELS))
  if (length(miss)) stop("resp_grid_wide: no ROW_LABELS entry for ", paste(miss, collapse = ", "))

  # facet_grid frees y BY ROW, so putting metrics on columns and asking for
  # free_y forces six different quantities onto one scale per region and
  # squashes every small-range metric. Each metric therefore gets its own plot,
  # with its own y scale, and patchwork lays them out side by side.
  one_col <- function(m, show_y_strip) {
    rg <- responsive_regions(m)
    d  <- plot_dat %>% transmute(Region, BIBI, BIBI_Class, y = .data[[m]])
    ggplot(d, aes(BIBI, y)) +
      { if (length(rg)) geom_smooth(data = d %>% filter(Region %in% rg),
          aes(colour = Region, fill = Region), method = "lm", formula = y ~ x,
          se = TRUE, linewidth = 0.5, alpha = 0.18, show.legend = FALSE) else NULL } +
      geom_point(aes(colour = Region, shape = BIBI_Class), size = 1.0, alpha = 0.9) +
      facet_grid(Region ~ ., switch = "y") +
      scale_colour_manual(values = REGION_COLS, limits = REGIONS, guide = "none") +
      scale_fill_manual(values = REGION_COLS, limits = REGIONS, guide = "none") +
      scale_shape_manual(values = CLASS_SHAPES, name = "Condition class", drop = FALSE) +
      scale_x_bibi(NULL) + scale_y_continuous(n.breaks = 3, labels = NEG_LAB) +
      labs(y = NULL, title = gsub("\n", " ", ROW_LABELS[[m]])) +
      theme_ch2 +
      theme(plot.title = element_text(size = FIG_BASE_PT - 2.5, face = "plain",
                                      hjust = 0.5, margin = margin(b = 2)),
            strip.placement = "outside",
            strip.text.y.left = if (show_y_strip)
              element_text(angle = 90, size = FIG_BASE_PT - 2) else element_blank(),
            strip.background.y = if (show_y_strip) element_rect(fill = "grey93", colour = NA)
                                 else element_blank(),
            panel.spacing = unit(3, "pt"),
            axis.text = element_text(size = FIG_BASE_PT - 3.5))
  }
  plots <- purrr::imap(metrics, function(m, k) one_col(m, k == 1))
  p <- patchwork::wrap_plots(plots, nrow = 1, guides = "collect") +
    patchwork::plot_annotation(caption = BIBI_AXIS_LABEL,
      theme = theme(plot.caption = element_text(hjust = 0.5, size = FIG_BASE_PT - 1,
                                                margin = margin(t = 2))))
  p & theme(legend.position = "bottom")
}

# One tolerance row, on relative abundance -- the weighting the index itself
# uses. Richness and Shannon cover the other two structural axes, and the two
# totals cover the currencies. The biomass-weighted HBI stays in the full screen (Table S4 part A).
if (WORKING_FIGS) {
  save_fig(resp_grid_wide(c("Richness_alpha", "Shannon", "HBI_rel_abund",
                            "Total_Density", "Total_Biomass")),
           "slide_structure_wide", 12.2, 4.35)
  save_fig(resp_grid_wide(c("FDiv_den", "FDis_den", "FR_RD_den",
                            "Gini_DB_DOM", "FR_DB_RD", "FR_RD_bio")),
           "slide_function_wide", 12.2, 4.35)
}

# The row lists are named because step 18 checks that FIGS3_METRICS exists (the scatter caption built from it is no
# longer written). Both grids are saved on every run although no manuscript carries them.
FIGS3_METRICS <- c("Richness_alpha", "Shannon", "HBI_rel_abund", "Total_Density", "Total_Biomass")
FIGS4_METRICS <- c("FDiv_den", "FDis_den", "FR_RD_den", "Gini_DB_DOM", "FR_DB_RD", "FR_RD_bio")
save_fig(resp_grid(FIGS3_METRICS), "figS3_structure_scatters", FIG_W_DOUBLE, 6.3)
save_fig(resp_grid(FIGS4_METRICS), "figS4_function_scatters", FIG_W_DOUBLE, 7.0)

say("STEP 06 complete")
