# =============================================================================
# 00_config.R — paths, constants, shared helpers
# =============================================================================
# WHY THIS FILE EXISTS
#   The thresholds, labels and input paths the pipeline uses are declared here
#   (the data folder itself in 00_paths.R). If an input moves, you change one line.
#
#   run_all.R is READ-ONLY on all canonical inputs. It never writes into the
#   data folder (only the data steps 00a to 00d do). Everything it produces goes to output/.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(purrr); library(stringr); library(vegan); library(ggplot2)
})

# ---- Project root ------------------------------------------------------------
# Resolved from this file's location, so the pipeline runs from any working
# directory. CH2_ROOT is the folder above the pipeline (used only in the project tree).
if (!exists("PIPE_DIR")) {
  PIPE_DIR <- tryCatch(
    dirname(dirname(normalizePath(sys.frame(1)$ofile))),
    error = function(e) normalizePath(".")
  )
}
# CH2_ROOT, DATA, RA, H1 and REPO_LAYOUT come from one file shared with the three data steps (00c, 00a, 00b).
# In the project tree DATA is <CH2_ROOT>/Data/Data; in the public repository it is <PIPE_DIR>/data (27 Sep 2026).
source(file.path(PIPE_DIR, "R", "00_paths.R"))

OUT       <- file.path(PIPE_DIR, "output")
OUT_TAB   <- file.path(OUT, "tables")
OUT_FIG   <- file.path(OUT, "figures")
OUT_LOG   <- file.path(OUT, "logs")
for (d in c(OUT, OUT_TAB, OUT_FIG, OUT_LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Canonical inputs --------------------------------------------------------
# Every one of these is post-Decapoda (regenerated 2026-08-06).
IN <- list(
  # Every screened metric, one row per site.
  # 18 Sept 2026: the redundancy family (RaoQ_gower_*, FR_RD_*, FR_DB_RD) was computed upstream on the BINARY
  # trait coding without group weights, while Supplement S1.5 describes the fuzzy group-weighted Gower for every
  # functional metric. R/00b_redundancy_geometry.R recomputes that family in the fuzzy geometry and writes
  # H1_site_summary_macro_metrics_fuzzyRD.csv (all other columns identical). Select it with CH2_RD_GEOMETRY=fuzzy.
  site_metrics = if (identical(Sys.getenv("CH2_RD_GEOMETRY"), "fuzzy"))
                   file.path(RA, "tables", "H1", "H1_site_summary_macro_metrics_fuzzyRD.csv")
                 else file.path(RA, "tables", "H1", "H1_site_summary_macro_metrics.csv"),

  # B-IBI. *** READ THE WARNING BELOW BEFORE CHANGING THIS PATH. ***
  # CH2_BIBI_MEDIAN=1 switches to the rarefied median score for every typology,
  # which is the basis of the published analysis. Without it the observed
  # chessie_bibi file is read, which is NOT the published basis; in the repository
  # run_all.R refuses to run without the flag.
  bibi         = if (nzchar(Sys.getenv("CH2_BIBI_MEDIAN")))
                   file.path(RA, "tables", "H1", "H1_BIBI_site_with_classes_MEDIAN.csv")
                 else file.path(RA, "tables", "H1", "H1_BIBI_site_with_classes.csv"),

  # Site x taxon matrices (sites by the trait-eligible taxa), exactly as passed to FD::dbFD.
  den_mat      = file.path(DATA, "den_mat_fd.rds"),
  bio_mat      = file.path(DATA, "bio_mat_fd.rds"),

  # Land-use gradient (PC1 of a within-region land-use PCA), one row per site.
  landuse      = file.path(RA, "env", "tables", "H2_gradient_scores_by_typology.csv"),

  # Taxon x trait-modality matrix (41 modalities) for steps 09, 09b, 11 and 15. Twardochleb et al. (2021) fuzzy traits, 22 Sep 2026:
  # the fuzzy file with genus melds and family gap fill; a group with no record anywhere in the family is NA.
  # Until 22 Sep 2026 this was traits_sel.rds, the legacy binary file, which is no longer a trait source.
  traits       = file.path(DATA, "traits_fuzzy_proxy_main.rds"),

  # --- FD provenance (step 13) -------------------------------------------
  # THESE, not traits_sel.rds, are the matrices the published FD metrics were
  # computed from. v1 ran with ANALYSIS_VERSION = "main", so pick_trait_file()
  # resolves to the *_proxy_main pair. Using traits_sel.rds instead reproduces
  # RaoQ at only r = 0.81. See 13_fd_provenance.R.
  traits_fuz   = file.path(DATA, "traits_fuzzy_proxy_main.rds"),
  traits_bin   = file.path(DATA, "traits_binary_proxy_main.rds"),

  # Block-balanced Gower weights (1 / block size), saved upstream with the trait matrices.
  trait_blocks = file.path(DATA, "trait_blocks.rds")
)

# Taxon tolerance values (step 09b). One value per taxon: the family-level Chessie B-IBI value the
# site HBI was computed from upstream, so the check uses exactly the numbers behind HBI. Classes follow
# the index's own metric definitions (Buchanan et al. 2011, Table 6): sensitive TV < 3, tolerant TV > 7.
# Added 21 September 2026.
IN$hbi_tv <- file.path(DATA, "H1_HBI_taxon_contributions.csv")
TV_SENSITIVE_BELOW <- 3
TV_TOLERANT_ABOVE  <- 7

# -----------------------------------------------------------------------------
# *** B-IBI SOURCE: read this before changing IN$bibi ***
#
# Every B-IBI value is the ICPRB release value: median_score in the three scoring files
# (Macros_Typ-1_..., Macros_Typ-2_..., Typ3_ChessieB-IBI_Scores_Metrics.csv in the data folder), and every class is the
# ICPRB median_rating. H1_BIBI_site_with_classes_MEDIAN.csv is written from them by
# build_bibi_site_classes.R, which rescales nothing and checks each class against its own bioregion's
# cut-points; step 01 stops unless every site equals its release value (tolerance 1e-6).
#
# Corrected 25 Sep 2026. An earlier version of this note called NORW's 0.3487 a factor-of-100 error and 34.87
# correct. It is the reverse: 0.3487 is NORW's score on the 0-100 scale (ICPRB mean over 100 rarefactions 1.81,
# sd 1.83, very poor in 100 of 100; it reproduces exactly from the MAC metric rules), and 34.87 came from
# an upstream per-value rescale, which multiplied any score <= 1.5 by 100. Never rescale a B-IBI value.
# -----------------------------------------------------------------------------

# ---- Analysis constants ------------------------------------------------------
REGIONS      <- c("Valley & Ridge", "Piedmont", "Coastal Plain")
RHO_MIN      <- 0.30    # screening: |rho| threshold
FDR_MAX      <- 0.10    # screening: BH-adjusted p threshold
P_EXPLORATORY<- 0.05    # exploratory tier: raw p threshold when FDR fails
N_PERM       <- 999     # PERMANOVA / BETADISPER permutations
N_PERM_ANOSIM<- 9999    # ANOSIM permutations (matches Methods 2.5)
SEED         <- 20260825
N_PERM_FC    <- 9999    # step 15 fourth-corner and RLQ permutations (22 Sep 2026). Until then step 15
                        # used N_PERM = 999, at which the states passing FDR 0.10 changed with the seed alone.
N_SEEDS_FC   <- 5       # the per-block fourth-corner (legacy Table 4 and Table S16, companion trait-state manuscript) runs under
                        # seeds SEED to SEED + 4; a state passes only when it meets FDR_MAX under every seed (22 Sep 2026).
                        # A state that passes under some seeds and not others is borderline, never passing.
N_PERM_SES   <- 999     # step 17 null model (trait labels shuffled); declared here since 27 Sep 2026, value unchanged
SES_RETAIN_SURVIVES <- 0.7   # step 17 verdict: share of the raw coefficient the SES coefficient keeps for "trait signal survives"
SES_RETAIN_PARTLY   <- 0.3   # and for "partly evenness-driven"; below it the signal is read as evenness, not traits

# -----------------------------------------------------------------------------
# METHOD CHOICES THAT DIFFER FROM v1 (the earlier version of this analysis) --
# this is what makes this a v2. v1 outputs are not part of this repository.
# -----------------------------------------------------------------------------

# TRANSFORM_PRIMARY -- square root, because that is what Methods 2.5 describes:
#   "Bray-Curtis dissimilarity on square-root-transformed density and biomass
#    matrices".
# v1's published within-region ANOSIM values were computed WITHOUT the
# transformation. v2 runs the analysis as described. Step 04 also runs the
# untransformed version and reports it side by side, so the difference between
# v1 and v2 is documented rather than silent.
#
# Consequence to be aware of: under the square-root transformation the Piedmont
# separates by condition class under BOTH weightings, not biomass only.
# See 04_multivariate_within_region.R.
TRANSFORM_PRIMARY <- "sqrt"
TRANSFORM_COMPARE <- "none"   # kept purely to reconcile against v1

# BH_SCOPE -- how wide the Benjamini-Hochberg family is.
#   "within_region"   : family = metrics in one family, in one region (v1 behavior,
#                       reproduces every published FDR value exactly)
#   "across_regions"  : family = metrics in one family, across all three regions
#                       (what Methods 2.5 literally says: "within each
#                       metric-hypothesis family ... across regions")
#
# These are NOT equivalent and the choice changes which metrics pass screening.
# Unlike the transformation, this one is a genuine judgment call rather than a
# discrepancy with a right answer. The paper uses within_region (v1 behavior);
# step 03 computes both and reports every metric whose tier would change.
BH_SCOPE <- "within_region"

# The screened metrics, grouped into their four BH families.
# BH correction is applied within Family x Region (see add_bh_within_family).
METRIC_FAMILIES <- list(
  alpha      = c("Richness_alpha","Shannon","HBI_rel_abund","HBI_bio"),
  structure  = c("Total_Density","Total_Biomass"),
  fd         = c("FRic_den","FEve_den","FEve_bio","FDiv_den","FDiv_bio",
                 "FDis_den","FDis_bio","RaoQ_den","RaoQ_bio"),
  redundancy = c("FR_RD_den","FR_RD_bio","FR_SF_den_log10","FR_DB_RD")
)

# NEW in v2: the concentration family (see 02_concentration_metrics.R).
# It is a SEPARATE BH family by design, so adding it cannot move any published
# FDR value in the four families above.
CONCENTRATION_FAMILY <- c("Gini_den","Gini_bio","N50_den","N50_bio","Gini_DB_DOM")

# ---- HBI is reported once, in its conventional form --------------------------
# HBI is a relative-abundance-weighted tolerance mean, and areal density is a
# constant multiple of relative abundance within a site, so HBI_den and
# HBI_rel_abund are the SAME statistic: identical rank correlations to 16
# significant figures in all three regions. Carrying both put a duplicate inside
# the BH denominator and made the structural correction conservative by one test,
# and it reported one metric as two. Relative abundance is the form Hilsenhoff
# (1988) defines, so that is the one screened.
#
# HBI_den is therefore removed everywhere -- family, dictionary, figures and
# verification -- rather than merely excluded from the screen.
#
# This CHANGES published FDR values in the alpha family. The only tier change is
# Coastal Plain Shannon, 0.1163 -> 0.0931, which now meets the 0.10 criterion.
# Decided 8 September 2026.

# -----------------------------------------------------------------------------
# PUBLICATION FIGURE AND TABLE SPECIFICATION
# -----------------------------------------------------------------------------
# Target: Freshwater Biology (Wiley). These are the numbers that matter, and
# getting them wrong is why the v1 figures did not survive scaling: a figure
# drawn 8.5 x 11 in and placed at 174 mm is reduced to 79% of its drawn size,
# which turns 10 pt axis text into 8 pt and 8 pt annotation into 6 pt. Draw at
# FINAL SIZE and the type stays the size you set.
#
#   single column   80 mm  = 3.15 in
#   double column  174 mm  = 6.85 in
#   max height     230 mm  = 9.05 in  (leaves room for the caption)
#   raster/line art 600 dpi
#   minimum type    7 pt at final size; body of figure at 8 pt
#   panel labels    (a) (b) (c), lower case, parentheses, bold, top left
FIG_W_SINGLE <- 3.15
FIG_W_DOUBLE <- 6.85
FIG_H_MAX    <- 9.05
FIG_DPI      <- 600
FIG_BASE_PT  <- 9.5   # serif needs more size than the sans default did

# Axis tick labels default to an ASCII hyphen for negatives, while the metric
# labels inside the same figures carry a true minus (U+2212) -- as does every
# coefficient in the manuscript. NEG_LAB makes the whole figure set agree.
# Added 10 September 2026 after the pre-submission figure audit.
NEG_LAB <- scales::label_number(style_negative = "minus", accuracy = NULL)

# Path to cairosvg, used to rebuild a .pdf from the .svg when R's own cairo
# device will not load. Empty string means "not installed", and save_fig then
# falls through to the lossy base device with a warning. Set CH2_CAIROSVG to
# override. See save_fig for why this route exists.
# (The authored figure captions live in the manuscript and are not read by the
# pipeline; until 27 Sep 2026 step 18 read them from a retired markdown export.)
CAIROSVG <- local({
  p <- Sys.getenv("CH2_CAIROSVG", "")
  if (nzchar(p) && file.exists(p)) return(p)
  for (cand in c(Sys.which("cairosvg")[[1]], "/opt/homebrew/bin/cairosvg",
                 "/usr/local/bin/cairosvg"))
    if (nzchar(cand) && file.exists(cand)) return(cand)
  ""
})

# Okabe-Ito, color-blind safe and distinguishable in grayscale. Regions keep a
# fixed color everywhere so a reader learns them once.
# ONE HUE FAMILY, ONE MEANING. Color encodes REGION and nothing else, anywhere
# in the figure set. Every other two-level split -- condition class, turnover vs
# nestedness -- is encoded by LUMINANCE (PAL_BINARY), which survives grayscale
# printing and every form of color-vision deficiency by construction, and which
# cannot be confused with a region.
#
# REGION HUE ASSIGNMENT IS FIXED BY THE STUDY-AREA MAP and must not be changed:
# Shenandoah / Valley & Ridge = GREEN, SE PA-MD Piedmont = BLUE,
# Delmarva / Coastal Plain = ORANGE, Urban/Developed = PURPLE (Chapter 3 only).
#
# The VALUES are the saturated Okabe-Ito equivalents of those hues, not the
# muted tones sampled from the Chapter 3 map. Same visual identity, materially
# better color-vision performance: worst-case separation between any two
# regions under deuteranopia or protanopia rises from 27 to 76. The muted map
# green and orange collapse under protanopia; these do not.
#
# PURPLE IS RESERVED for the Chapter 3 Urban/Developed typology and must never
# appear in these figures.
PAL_REGION <- c("Valley & Ridge" = "#009E73",   # map hue: green
                "Piedmont"       = "#0072B2",   # map hue: blue
                "Coastal Plain"  = "#D55E00")   # map hue: orange
PAL_URBAN  <- "#CC79A7"                          # reserved, Chapter 3 only

# Two-level ordinal splits. Dark = the first / higher-magnitude level.
PAL_BINARY <- c(dark = "#2B2B2B", light = "#B4B4B4")

PAL_CONDITION <- c("Low condition (bottom third)"  = unname(PAL_BINARY["light"]),
                   "High condition (top third)"    = unname(PAL_BINARY["dark"]))
PAL_COMPONENT <- c("Nestedness" = unname(PAL_BINARY["dark"]),
                   "Turnover"   = unname(PAL_BINARY["light"]))
PAL_CURRENCY  <- c("Density" = "#332288", "Biomass" = "#88CCEE")

# Condition class is ORDINAL and shared across steps. Shape carries it in the
# response figures; a sequential gray ramp carries it where shape cannot.
# The five ICPRB narrative ratings, taken from median_rating rather than binned
# by us. Ordered worst to best. Under CH2_BIBI_MEDIAN the input carries these;
# the four straddle bins remain for the default (chessie) run.
CLASS_LEVELS <- if (nzchar(Sys.getenv("CH2_BIBI_MEDIAN")))
  c("Very poor", "Poor", "Fair", "Good", "Excellent") else
  c("VeryPoor–Poor", "Poor–Fair", "Fair–Good", "Good–Excellent")
# REGION_COLS is an alias of PAL_REGION for the steps that use that name; region
# color is the ONLY use of hue in the figure set. Step 06 once defined its own
# Dark2 palette, which collided with the condition-class and beta-component colors.
REGION_COLS  <- PAL_REGION
CLASS_SHAPES <- if (nzchar(Sys.getenv("CH2_BIBI_MEDIAN")))
  c("Very poor" = 16, "Poor" = 17, "Fair" = 15, "Good" = 18, "Excellent" = 8) else
  c("VeryPoor–Poor" = 16, "Poor–Fair" = 17, "Fair–Good" = 15, "Good–Excellent" = 8)


#' The one figure theme. SERIF, to match the manuscript body face -- a
#' deliberate choice, overriding the earlier sans-serif default. The reason that default
#' existed was legibility: serif at 7-8 pt on a reduced figure is the first
#' thing to become unreadable. The base size is therefore raised to 9.5 pt to
#' pay for the switch, and every figure is drawn at final placed size so what
#' is set here is what prints.
FIG_FAMILY <- local({
  f <- try(systemfonts::system_fonts(), silent = TRUE)
  if (!inherits(f, "try-error") && any(f$family == "Times New Roman")) "Times New Roman" else "serif"
})
# The PNG and SVG devices reach the real system font through systemfonts. The
# PDF device does not: it carries its own PostScript font database, which knows
# the fourteen base fonts and nothing else, so "Times New Roman" is not found
# and any bold or italic run then fails outright with "invalid font type".
# Register the name against base PostScript Times, whose metrics Times New Roman
# was drawn to match. cairo_pdf would sidestep this, but cairo can fail to load
# on macOS builds without XQuartz (see save_fig) -- for every family, not just this one.
if (FIG_FAMILY == "Times New Roman" &&
    !"Times New Roman" %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts("Times New Roman" = grDevices::pdfFonts()$Times)
  grDevices::postscriptFonts("Times New Roman" = grDevices::postscriptFonts()$Times)
}
# Serif is wider than the sans face these panels were first drawn in, so the
# longer framework titles ran off the rightmost column. Wrap to a fixed TWO
# lines -- short titles get a padded blank second line -- so the full wording
# survives and every panel in a row still starts at the same height.
wrap_title <- function(x, width = 30) {
  w <- paste(strwrap(x, width = width), collapse = "\n")
  if (!grepl("\n", w)) w <- paste0(w, "\n")
  w
}

ggplot2::theme_update(text = ggplot2::element_text(family = FIG_FAMILY))

theme_pub <- function(base = FIG_BASE_PT) {
  ggplot2::theme_bw(base_size = base, base_family = FIG_FAMILY) +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      panel.grid.major  = ggplot2::element_line(linewidth = 0.25, colour = "grey88"),
      panel.border      = ggplot2::element_rect(linewidth = 0.4, colour = "grey35"),
      axis.text         = ggplot2::element_text(size = base - 1, colour = "grey20"),
      axis.title        = ggplot2::element_text(size = base),
      strip.background  = ggplot2::element_rect(fill = "grey94", colour = NA),
      strip.text        = ggplot2::element_text(size = base, face = "bold",
                                                margin = ggplot2::margin(2.5,2.5,2.5,2.5)),
      legend.key.size   = grid::unit(0.24, "cm"),
      legend.box.margin = ggplot2::margin(-4,0,0,0),
      legend.spacing.y  = grid::unit(0.02, "cm"),
      legend.box.spacing= grid::unit(0.12, "cm"),
      legend.text       = ggplot2::element_text(size = base - 1),
      legend.title      = ggplot2::element_text(size = base - 1),
      legend.margin     = ggplot2::margin(0,0,0,0),
      plot.tag          = ggplot2::element_text(size = base + 1, face = "bold"),
      plot.title        = ggplot2::element_text(size = base, face = "bold"),
      plot.caption      = ggplot2::element_text(size = base - 2, hjust = 0,
                                                colour = "grey30"),
      panel.spacing     = grid::unit(0.30, "lines"),
      plot.margin       = ggplot2::margin(1,2,1,1))
}

#' The B-IBI axis, bounded to the full index range.
#'
#' The index runs 0-100 and a reader needs to see where the observed data sit
#' within that range. Letting ggplot auto-scale to the observed 6.2-88.6 hides
#' the bound and makes the three regions look as though they span different
#' scales when they do not. Breaks every 25 give the quartile grid without
#' crowding at 8 pt.
scale_x_bibi <- function(title = NULL) {
  list(
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                                expand = ggplot2::expansion(mult = 0.02)),
    ggplot2::labs(x = title))
}

#' Panel labels the journal way: (a), (b), (c), bold, top left.
# Do NOT set plot.tag.position to npc coordinates here. With a panel title
# present the tag is drawn into the title strip and is clipped -- the first
# render of Figure 4 lost both labels that way. patchwork's default position
# places the tag above and left of the title, which is what the journal wants.
# Takes the assembled patchwork and returns it tagged. It must work this way:
# `&` applies a theme across every subplot of a PATCHWORK, so it has to be
# applied to the whole assembled object. A helper that returns
# `plot_annotation(...) & theme(...)` and is then added with `+` silently
# produces no tags at all -- that is how Figure 4 lost its labels twice.
tag_panels <- function(p) {
  p + patchwork::plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
    ggplot2::theme(plot.tag = ggplot2::element_text(size = FIG_BASE_PT + 2,
                                                    family = FIG_FAMILY,
                                                    face = "bold", hjust = 0, vjust = 1))
}

#' Round every numeric column for a manuscript-facing table. Raw CSVs carry
#' 17 significant digits, which is not a table anyone can read or paste.
fmt_table <- function(df, digits = 3) {
  dplyr::mutate(df, dplyr::across(dplyr::where(is.numeric), ~ round(.x, digits)))
}

# -----------------------------------------------------------------------------
# THE METRIC DICTIONARY -- ONE SOURCE OF TRUTH FOR EVERY NAME IN THE PAPER
# -----------------------------------------------------------------------------
# Added 2026-08-27 after figure labels drifted away from the manuscript's own
# vocabulary. Three things had gone wrong and all three were avoidable:
#   - panel titles name the family as the field names it (25 Sep 2026:
#     functional diversity / redundancy / inequality; totals; Methods 2.4);
#   - panels said "numerical" and "energetic" where the standing rule is one
#     word per concept: DENSITY and BIOMASS;
#   - the x axis read "B-IBI score (0-100 index...)" rather than naming the
#     gradient it represents.
#
# Every figure label, table row and caption now comes from this table, so a
# name can only be wrong here, once, rather than in six places. Step 18 also
# emits it as a manuscript table so the reader sees the same definitions.
#
# Columns: key       pipeline column name
#          property  functional diversity / redundancy / inequality / structure / totals
#          name      the established name, as used in Methods 2.4
#          symbol    the abbreviation the manuscript uses
#          currency  Density / Biomass / Relative abundance (HBI only) / invariant / difference   (the weighting; 'currency' is the retired word, kept as a column name only)
#          axis      y-axis label, with units where they exist
#          equation  where one is needed to be unambiguous
#          citation  the source for the index, narrowed to what the paper supports (checked against
#                    the full texts, 28 Sep 2026; the Source column of Table S2 part C, legacy Table S7)
METRIC_DICT <- tibble::tribble(
  ~key, ~property, ~name, ~short, ~meaning, ~symbol, ~currency, ~axis, ~equation, ~citation, ~y_acc,
  "Richness_alpha", "structure", "Taxonomic richness", "taxonomic richness", "Number of taxa present, each identified to the finest rank practicable; independent of weighting.", "S", "invariant", "S (taxa)", "", "count of taxa; no primary source", "10",
  "Shannon", "structure", "Shannon diversity", "Shannon diversity", "Evenness-weighted diversity of the community.", "H'", "Density", "H'", "H' = -sum(p_i ln p_i)", "Shannon (1948)", "0.5",
  "HBI_rel_abund", "structure", "Hilsenhoff Biotic Index", "Hilsenhoff Biotic Index", "Weighted mean pollution tolerance; lower values indicate less tolerant taxa. The conventional relative-abundance form; weighting by areal density reproduces it exactly.", "HBI", "Relative abundance", "HBI", "", "Hilsenhoff (1988)", "0.5",
  "HBI_bio", "structure", "Hilsenhoff Biotic Index", "Hilsenhoff Biotic Index", "Weighted mean pollution tolerance; lower values indicate less tolerant taxa.", "HBI", "Biomass", "HBI", "", "Hilsenhoff (1988)", "0.5",
  "Total_Density", "totals", "Total density", "total density", "Areal abundance of all quantitatively sampled taxa.", "-", "Density", "Density (ind. m\u207b\u00b2)", "", "Moulton et al. (2000)", "0",
  "Total_Biomass", "totals", "Total biomass", "total biomass", "Areal dry mass of all quantitatively sampled taxa.", "-", "Biomass", "Biomass (mg DM m\u207b\u00b2)", "", "Benke et al. (1999)", "0",
  "FRic_den", "functional diversity", "Functional richness", "functional richness", "Volume of trait space enclosed by the taxa present; presence-based, so the same with both weightings.", "FRic", "invariant", "FRic", "convex hull volume", "Villeger et al. (2008)", "0.05",
  "FEve_den", "functional diversity", "Functional evenness", "functional evenness", "Regularity with which weighted taxa are spaced through the occupied trait space.", "FEve", "Density", "FEve", "regularity along the MST", "Villeger et al. (2008)", "0.05",
  "FEve_bio", "functional diversity", "Functional evenness", "functional evenness", "Regularity with which weighted taxa are spaced through the occupied trait space.", "FEve", "Biomass", "FEve", "regularity along the MST", "Villeger et al. (2008)", "0.05",
  "FDiv_den", "functional diversity", "Functional divergence", "functional divergence", "How peripheral the most abundant taxa are within the occupied trait space.", "FDiv", "Density", "FDiv", "deviation from mean distance to hull center", "Villeger et al. (2008)", "0.05",
  "FDiv_bio", "functional diversity", "Functional divergence", "functional divergence", "How peripheral the most abundant taxa are within the occupied trait space.", "FDiv", "Biomass", "FDiv", "deviation from mean distance to hull center", "Villeger et al. (2008)", "0.05",
  "RaoQ_den", "functional diversity", "Rao's quadratic entropy", "Rao's quadratic entropy", "Weighted mean pairwise trait dissimilarity among co-occurring taxa.", "RaoQ", "Density", "RaoQ", "Q = sum_ij p_i p_j d_ij", "Rao (1982); Botta-Dukat (2005)", "0.01",
  "RaoQ_bio", "functional diversity", "Rao's quadratic entropy", "Rao's quadratic entropy", "Weighted mean pairwise trait dissimilarity among co-occurring taxa.", "RaoQ", "Biomass", "RaoQ", "Q = sum_ij p_i p_j d_ij", "Rao (1982); Botta-Dukat (2005)", "0.01",
  "FDis_den", "functional diversity", "Functional dispersion", "functional dispersion", "Weighted mean distance of taxa to the weighted community centroid.", "FDis", "Density", "FDis", "mean distance to weighted centroid", "Laliberte & Legendre (2010)", "0.01",
  "FDis_bio", "functional diversity", "Functional dispersion", "functional dispersion", "Weighted mean distance of taxa to the weighted community centroid.", "FDis", "Biomass", "FDis", "mean distance to weighted centroid", "Laliberte & Legendre (2010)", "0.01",
  "FR_RD_den", "redundancy", "Overlap redundancy", "overlap redundancy", "Fraction of Simpson diversity not expressed as functional dissimilarity; the degree to which taxa share roles.", "FR_RD", "Density", "FR_RD", "FR_RD = 1 - Q/D", "Ricotta et al. (2016), eqn 4, who name it functional redundancy; after de Bello et al. (2007), whose unnormalized form is D - Q", "0.05",
  "FR_RD_bio", "redundancy", "Overlap redundancy", "overlap redundancy", "Fraction of Simpson diversity not expressed as functional dissimilarity; the degree to which taxa share roles.", "FR_RD", "Biomass", "FR_RD", "FR_RD = 1 - Q/D", "Ricotta et al. (2016), eqn 4, who name it functional redundancy; after de Bello et al. (2007), whose unnormalized form is D - Q", "0.02",
  "FR_SF_den_log10", "redundancy", "Space-filling redundancy", "space-filling redundancy", "Packing density of trait-eligible taxa within the occupied hull.", "FR_SF", "invariant", "log10 FR_SF", "FR_SF = S_taxa / FRic", "adapted from Mouillot et al. (2014), with FRic (Villeger et al., 2008) in place of the number of functional entities", "0.5",
  "FR_DB_RD", "redundancy", "Redundancy divergence", "redundancy divergence", "Difference between biomass- and density-weighted overlap redundancy; positive values indicate greater trait sharing among the taxa that make up most of the biomass.", "FR_DB_RD", "difference", "FR_RD (biomass \u2212 density)", "FR_RD(biomass) - FR_RD(density)", "this study (extends Cummins et al., 2022)", "0.05",
  "Gini_den", "inequality", "Gini coefficient", "Gini coefficient", "Inequality of the taxon share vector; higher values indicate individuals concentrated in fewer taxa.", "Gini", "Density", "Gini", "", "Gini (1912); applied to the sizes of plant individuals by Weiner & Solbrig (1984), here to taxon shares", "0.05",
  "Gini_bio", "inequality", "Gini coefficient", "Gini coefficient", "Inequality of the taxon share vector; higher values indicate biomass concentrated in fewer taxa.", "Gini", "Biomass", "Gini", "", "Gini (1912); applied to the sizes of plant individuals by Weiner & Solbrig (1984), here to taxon shares", "0.05",
  "N50_den", "inequality", "Taxa making up half the total", "taxa making up half the total", "Number of taxa needed to reach half the community total; higher values indicate a more even spread among taxa.", "N50", "Density", "N50 (taxa)", "", "constructed here", "5",
  "N50_bio", "inequality", "Taxa making up half the total", "taxa making up half the total", "Number of taxa needed to reach half the community total; higher values indicate a more even spread among taxa.", "N50", "Biomass", "N50 (taxa)", "", "constructed here", "5",
  "Gini_DB_DOM", "inequality", "Concentration divergence", "concentration divergence", "Difference between biomass- and density-weighted Gini coefficients; positive values indicate biomass spread more unequally among taxa than density.", "Gini_DB_DOM", "difference", "Gini (biomass \u2212 density)", "Gini(biomass) - Gini(density)", "this study", "0.05"
)

#' Look up one field for one metric. Fails loudly rather than silently
#' returning the raw column name, which is how the wrong labels survived.
mdict <- function(key, field) {
  r <- METRIC_DICT[METRIC_DICT$key == key, ]
  if (!nrow(r)) stop("mdict(): '", key, "' is not in METRIC_DICT. Add it there, ",
                     "do not hard-code a label at the call site.")
  r[[field]][1]
}

#' Panel title: the established name, prefixed by the framework property.
#' "Redundancy: overlap redundancy (density)", never "Redundancy: numerical".
mdict_title <- function(key) {
  p <- mdict(key, "property"); cu <- mdict(key, "currency")
  # `short` is stored with the exact casing wanted after the colon. Do NOT
  # tolower() the name here -- that produced "Geometry: rao's quadratic entropy".
  ttl <- if (p %in% c("functional diversity", "redundancy", "inequality"))   # 25 Sep 2026: the families as the field names them
           sprintf("%s: %s", paste0(toupper(substring(p,1,1)), substring(p,2)),
                   mdict(key, "short"))
         else mdict(key, "name")
  # Standing-stock names already carry the currency ("Total biomass"), so do
  # not append it again.
  if (cu %in% c("Density", "Biomass") && p != "totals")
    sprintf("%s (%s)", ttl, tolower(cu)) else ttl
}

#' Y-axis label: symbol plus currency, with units where they exist.
mdict_axis <- function(key) {
  a <- mdict(key, "axis"); cu <- mdict(key, "currency")
  if (cu %in% c("Density", "Biomass") && !grepl("m\u207b\u00b2", a))
    sprintf("%s (%s)", a, tolower(cu)) else a
}

#' Y scale with the same look in every facet of a metric.
#'
#' Free y ranges are right -- a shared range let one Valley & Ridge outlier
#' flatten the other two regions. But free ranges also let ggplot pick its own
#' breaks per facet, which gave HBI as 3/4/5/6 in Valley & Ridge, 4.0/4.5/5.0
#' in the Piedmont and 5/6 in the Coastal Plain: one metric, three formats.
#' This snaps every facet to the same rounding step and the same number of
#' breaks, so the values differ but the formatting does not.
# -----------------------------------------------------------------------------
# FIXED Y DOMAINS -- one axis per metric, shared across regions AND currencies
# -----------------------------------------------------------------------------
# Added 2026-08-27 after the axes were still being read off each panel's own
# extremes. Labels like 0.72 / 0.81 / 0.90 / 0.99 are visibly derived from the
# data rather than chosen, they differ between the two currencies of one
# metric, and they invite the obvious review question: why does the axis not
# start at zero, and where does this sit in the metric's possible range?
#
# Each metric now declares a round domain that brackets its observed range
# across every panel it appears in, with a round step. Consequences:
#   - FR_RD density and FR_RD biomass share one axis, so the two backup panels
#     are directly comparable instead of merely similar-looking;
#   - metrics with a meaningful zero START at zero (RaoQ, FDis, FEve, FRic,
#     richness, N50, HBI), which is the honest answer to the review question;
#   - bounded indices whose data occupy a narrow high band (FDiv, FR_RD, Gini)
#     start at a round value below the data rather than at zero, because a
#     0-1 axis compresses the within-region trend to invisibility. The [0,1]
#     bound is stated in the caption instead.
#
# Each domain is c(lo, hi, step).
# Domains are DERIVED from the data, not hand-listed, so adding sites or a
# metric cannot silently leave an axis stale or clip a point.
#
# THE RULE
#   non-negative metric : axis runs 0 to a round value at or above the maximum
#   signed metric       : axis runs round(min) to round(max) and always
#                         includes 0 as a labeled break
# Either way the axis shows the full scale the data sit within -- the same
# thing the B-IBI axis does by running 0 to 100 while the data span 6 to 89 --
# so a reader can see where a region falls in the possible range rather than
# only how it varies internally.
.nice_ceiling <- function(m) {
  if (!is.finite(m) || m == 0) return(1)
  k <- 10^floor(log10(abs(m)))
  c0 <- c(1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10)
  cand <- c0 * k
  cand[which(cand >= abs(m) - 1e-12)[1]] * sign(m)
}

#' One domain per metric symbol: c(lo, hi, step).
compute_y_domains <- function(data) {
  out <- list()
  for (sym in setdiff(unique(METRIC_DICT$symbol), "-")) {
    ks <- intersect(METRIC_DICT$key[METRIC_DICT$symbol == sym], names(data))
    if (!length(ks)) next
    v <- unlist(data[ks]); v <- v[is.finite(v)]
    if (length(v) < 2) next
    lo0 <- min(v); hi0 <- max(v)
    mag  <- 10^floor(log10(max(abs(c(lo0, hi0)))))
    steps <- sort(c(1, 2, 2.5, 5, 10) * rep(mag * c(0.01, 0.1, 1), each = 5))
    pick <- NULL
    for (st in steps) {
      # anchoring both ends to multiples of the step guarantees a zero label
      lo <- if (lo0 >= 0) 0 else -ceiling((-lo0)/st - 1e-12) * st
      hi <- ceiling(hi0/st - 1e-12) * st
      nb <- round((hi - lo)/st) + 1
      if (nb >= 4 && nb <= 6) { pick <- c(lo, hi, st); break }
    }
    if (is.null(pick)) {                      # nothing in range: fall back
      hi <- .nice_ceiling(hi0)
      lo <- if (lo0 >= 0) 0 else -.nice_ceiling(-lo0)
      pick <- c(lo, hi, (hi - lo)/4)
    }
    out[[sym]] <- pick
  }
  # standing stock carries no symbol of its own
  for (k in c("Total_Density", "Total_Biomass")) {
    if (!k %in% names(data)) next
    v <- data[[k]]; v <- v[is.finite(v)]
    hi <- .nice_ceiling(max(v))
    out[[if (k == "Total_Density") "density" else "biomass"]] <- c(0, hi, hi/4)
  }
  out
}

# Populated by step 06 from the assembled plotting data.
Y_DOMAIN <- list()

#' Exactly N y-axis breaks, uniformly spaced, on the label grid.
#'
#' WHY THE COUNT RULE IS HAND-WRITTEN RATHER THAN A PACKAGE CALL. Checked:
#'   scales::breaks_log(n = 4)       890-9700  -> 1,000 3,000          (n = 2)
#'                                   235-3390  -> 300 1,000 3,000      (n = 3)
#'   scales::breaks_extended(n = 4)  3.05-6.60 -> 4 5 6                (n = 3)
#'                                   4.40-6.60 -> 4.5 5.0 5.5 6.0 6.5  (n = 5)
#' In every standard implementation -- ggplot/labeling (Wilkinson extended),
#' base pretty(), matplotlib MaxNLocator -- n is a TARGET, not a guarantee,
#' because a fixed count and conventionally round values are in tension. No
#' library ships a fixed-count breaker, so this rule has to live here.
#' Formatting and log ladders, which packages DO solve well, are delegated to
#' scales rather than hand-rolled.
#'
#' THE RULE, after three earlier attempts got half of it:
#'   - the NUMBER of labels is fixed (default 4) in every panel;
#'   - the DECIMALS are fixed per metric, so a metric reads the same wherever
#'     it appears, including its other currency;
#'   - the STEP therefore varies between panels, because the ranges genuinely
#'     differ. Holding step and count both constant is impossible with free
#'     ranges, and count is the one a reader notices.
#'
#' Earlier failures worth not repeating: thinning a fine grid gave uneven
#' spacing (30/40/60/70); snapping the step upward gave 2, 3 and 2 labels
#' across one metric; and a step of 0.125 printed at two decimals rendered as
#' -0.12/0.00/0.12/0.25, which looks irregular even though it is not.
#'
#' @param r     panel limits
#' @param grid  smallest representable label increment (10^-dec)
#' @param n     number of breaks to show
breaks_exact_n <- function(r, grid, n = 4) {
  span <- diff(r)
  if (!is.finite(span) || span <= 0 || !is.finite(grid) || grid <= 0)
    return(numeric(0))
  # THE RULE: n labels, evenly spaced, spanning as much of the axis as
  # possible -- what the B-IBI axis does with 0/25/50/75/100. Coverage beats
  # roundness; optimizing for round steps gave 0.60/0.70/0.80/0.90 in one
  # backup panel and 0.64/0.68/0.72/0.76 in the next.
  #
  # CRITICAL: ggplot passes the EXPANDED limits, and any break outside them is
  # silently dropped. Breaks must therefore be fitted INSIDE `r`, never
  # generated from a rounded-outward version of it -- doing that left only two
  # labels visible per panel.
  lo <- ceiling((r[1] - 1e-12)/grid) * grid
  hi <- floor((r[2] + 1e-12)/grid) * grid
  usable <- hi - lo
  if (!is.finite(usable) || usable <= 0) return(numeric(0))

  fits <- function(st) {
    if ((n - 1) * st > usable + 1e-12) return(NULL)
    a <- ceiling((lo - 1e-12)/st) * st                  # align to the step
    if (a + (n - 1) * st > hi + 1e-12)
      a <- lo + floor(((usable - (n - 1) * st)/2)/grid) * grid   # else center
    if (a < lo - 1e-12 || a + (n - 1) * st > hi + 1e-12) return(NULL)
    a + st * (0:(n - 1))
  }
  # Prefer a conventional step, but only one that still fills most of the axis.
  mult <- c(1, 1.5, 2, 2.5, 3, 4, 5, 6, 10, 15, 20, 25, 30, 40, 50,
            60, 100, 150, 200, 250, 300, 500, 1000)
  cand <- sort(mult * grid, decreasing = TRUE)
  cand <- cand[abs(cand/grid - round(cand/grid)) < 1e-9]
  for (st in cand) {
    b <- fits(st)
    if (!is.null(b) && (n - 1) * st >= 0.70 * usable) return(b)
  }
  # Otherwise the largest step on the grid that fits -- maximizes coverage.
  st <- floor((usable/(n - 1))/grid) * grid
  if (st < grid) st <- grid
  b <- fits(st)
  if (is.null(b)) seq(lo, hi, length.out = n) else b
}

#' Limits matched to the breaks, so the axis really does run bottom to top.
limits_for_breaks <- function(b, r) {
  if (!length(b)) return(NULL)
  c(min(b[1], r[1]), max(b[length(b)], r[2]))
}

#' Y scale for one metric: fixed break count, decimals fixed per metric.
scale_y_metric <- function(key, steps = NULL, n = 4) {
  # The derived domain is checked FIRST. Standing stock carries y_acc = "0",
  # so testing y_acc before the domain lookup returned a default scale and the
  # density and biomass panels silently lost their shared 0-to-maximum axis.
  sym <- mdict(key, "symbol")
  if (identical(sym, "-"))                       # standing stock has no symbol
    sym <- if (key == "Total_Density") "density" else
           if (key == "Total_Biomass") "biomass" else sym
  if (sym %in% names(Y_DOMAIN)) {
    d  <- Y_DOMAIN[[sym]]
    br <- seq(d[1], d[2], by = d[3])
    dd <- max(0, -floor(log10(d[3])))
    return(ggplot2::scale_y_continuous(
      limits = c(d[1], d[2]), breaks = br,
      expand = ggplot2::expansion(mult = 0.02),
      # big.mark must be given explicitly; the default in this locale is a
      # space, which printed "10 000" rather than "10,000".
      labels = scales::label_number(accuracy = 10^(-dd), big.mark = ",")))
  }

  acc <- suppressWarnings(as.numeric(mdict(key, "y_acc")))
  if (is.na(acc) || acc <= 0) return(ggplot2::scale_y_continuous(n.breaks = n))
  dec  <- max(0, -floor(log10(acc)))
  grid <- 10^(-dec)
  ggplot2::scale_y_continuous(
    breaks = function(limits) breaks_exact_n(range(limits, na.rm = TRUE), grid, n),
    # a 2% expansion, so the outermost labels sit close to the ends of the axis,
    # as 0 and 100 do on the B-IBI axis
    expand = ggplot2::expansion(mult = 0.02),
    # scales::label_number is the idiomatic formatter and will not print
    # float noise the way formatC does.
    labels = scales::label_number(accuracy = grid, big.mark = ","))
}

#' The same rule on a log10 axis, for the standing-stock panels.
#'
#' A log axis needs a CONVENTIONAL ladder. Interpolating evenly in log space
#' and rounding to two significant figures produced 839 / 1,899 / 4,499 /
#' 10,000, which are not axis numbers. Walk progressively finer 1-2-3-5 style
#' ladders and take the first that supplies enough labels inside the range.
scale_y_log_n <- function(n = 4) {
  ggplot2::scale_y_log10(
    breaks = function(limits) {
      r <- range(limits, na.rm = TRUE)
      if (!all(is.finite(r)) || r[1] <= 0 || diff(r) <= 0) return(numeric(0))
      ks      <- floor(log10(r[1])):ceiling(log10(r[2]))
      ladders <- list(c(1), c(1, 3), c(1, 2, 5), c(1, 2, 3, 5),
                      c(1, 1.5, 2, 3, 5, 7))
      pick <- NULL
      for (L in ladders) {
        # 10^k is computed as exp(k*log(10)) and lands a hair under the exact
        # power, so 5 * 10^3 becomes 4999.9999... and printed as "4,999".
        # Round the ladder to significant figures before it is ever formatted.
        v <- sort(signif(as.vector(outer(L, 10^ks)), 6))
        v <- v[v >= r[1] & v <= r[2]]
        if (length(v) >= n) { pick <- v; break }
        if (is.null(pick) || length(v) > length(pick)) pick <- v
      }
      if (!length(pick)) return(numeric(0))
      if (length(pick) > n) {
        drop  <- length(pick) - n
        front <- floor(drop/2)
        pick  <- pick[(front + 1):(front + n)]
      }
      pick
    },
    labels = scales::label_number(accuracy = 1, big.mark = ","))
}

# The gradient, named once. Every x axis uses this.
BIBI_AXIS_LABEL <- "Biological condition gradient (B-IBI)"

# ---- Shared helpers ----------------------------------------------------------

#' Map any of the region spellings used across the project to one canonical set.
as_region <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(
    x %in% c("1","VR","ValleyRidge","Valley & Ridge")        ~ "Valley & Ridge",
    x %in% c("3","PD","Piedmont","NorthernPiedmont")         ~ "Piedmont",
    x %in% c("2","CP","CoastalPlain","Coastal Plain")        ~ "Coastal Plain",
    TRUE ~ NA_character_
  )
  if (any(is.na(out))) stop("as_region(): unmapped region label(s): ",
                            paste(unique(x[is.na(out)]), collapse = ", "))
  factor(out, levels = REGIONS)
}

std_site <- function(x) toupper(trimws(as.character(x)))

#' Spearman rho + raw p for one metric in one region.
spearman_one <- function(bibi, metric, region_label, metric_label) {
  ok <- is.finite(bibi) & is.finite(metric)
  if (sum(ok) < 4) {
    return(tibble(Metric = metric_label, Region = region_label, n = sum(ok),
                  rho = NA_real_, p_raw = NA_real_))
  }
  ct <- suppressWarnings(cor.test(bibi[ok], metric[ok], method = "spearman",
                                  exact = FALSE))
  tibble(Metric = metric_label, Region = region_label, n = sum(ok),
         rho = unname(ct$estimate), p_raw = ct$p.value)
}

#' Benjamini-Hochberg, applied within each family WITHIN each region.
#'
#' Methods 2.5 reads "within each metric-hypothesis family ... across regions",
#' which is ambiguous. The canonical published FDR values settle it: the
#' correction is applied within Family x Region, with n = the number of metrics
#' in that family. Verified against every published FDR value, e.g.
#'   fd family, Piedmont, n = 9 : FDiv_den raw 0.0055 -> 0.0497  (published)
#'   structure,  Piedmont, n = 2 : Total_Biomass 0.0141 -> 0.0283 (published)
#'   redundancy, Coastal Plain, n = 4 : FR_RD_den 0.0005 -> 0.0020 (published)
#' Grouping by Family alone gives different numbers.
add_bh_within_family <- function(df, scope = BH_SCOPE) {
  if (scope == "within_region") {
    df %>% group_by(Family, Region) %>%
      mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>% ungroup()
  } else if (scope == "across_regions") {
    df %>% group_by(Family) %>%
      mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>% ungroup()
  } else stop("BH_SCOPE must be 'within_region' or 'across_regions'")
}

#' The screening rule from Methods 2.5.
add_tier <- function(df) {
  df %>% mutate(
    Tier = case_when(
      is.na(rho)                                      ~ "not testable",
      abs(rho) >= RHO_MIN & p_FDR <= FDR_MAX          ~ "strong",
      abs(rho) >= RHO_MIN & p_raw  <= P_EXPLORATORY   ~ "exploratory",
      TRUE                                            ~ "not responsive"
    )
  )
}

#' Save a figure as PNG + PDF, robustly.
#'
#' Device choice matters here and is not cosmetic. The default pdf() device
#' substitutes plain hyphens for en-dashes and minus signs, which mangles the
#' condition-class labels ("VeryPoor-Poor") and the "biomass - density" axis
#' titles. cairo_pdf handles UTF-8 -- but capabilities("cairo") can report TRUE
#' on macOS while the cairo DLL still fails to load because the X11 libraries
#' are missing. So: try cairo, fall back to pdf(), and say which was used.
save_fig <- function(p, name, w, h, dpi = FIG_DPI) {
  # Guard the journal limits. Drawing wider than a page means the figure is
  # scaled down at typesetting and every type size shrinks with it.
  if (w > FIG_W_DOUBLE + 0.01)
    say("  WARNING: %s is %.2f in wide; %s max is %.2f in (174 mm); type will shrink when it is scaled to fit",
        name, w, "Freshwater Biology", FIG_W_DOUBLE)
  if (h > FIG_H_MAX + 0.01)
    say("  WARNING: %s is %.2f in tall; max with caption is %.2f in (230 mm)",
        name, h, FIG_H_MAX)
  ggplot2::ggsave(file.path(OUT_FIG, paste0(name, ".png")), p,
                  width = w, height = h, dpi = dpi, device = ragg::agg_png)

  # SVG IS WRITTEN BEFORE THE PDF, and that order is load-bearing.
  # svglite goes through systemfonts, so U+2212 MINUS and U+2013 EN DASH
  # survive as themselves and text stays text rather than outlines. That makes
  # the SVG the faithful vector record -- and, when cairo is unavailable, the
  # thing the PDF is rebuilt from. It is also what lets a figure be
  # hand-adjusted for placement in Illustrator or Inkscape.
  svg_path <- file.path(OUT_FIG, paste0(name, ".svg"))
  try(suppressWarnings(ggplot2::ggsave(svg_path, p, width = w, height = h,
                       device = svglite::svglite)), silent = TRUE)

  # cairo failing is a WARNING, not an error, and ggsave then returns having
  # written nothing. So test the file, do not trust the call.
  pdf_path <- file.path(OUT_FIG, paste0(name, ".pdf"))
  if (file.exists(pdf_path)) unlink(pdf_path)
  try(suppressWarnings(ggplot2::ggsave(pdf_path, p, width = w, height = h,
                                       device = grDevices::cairo_pdf)), silent = TRUE)

  # FALLBACK 1 -- rebuild the PDF from the SVG.
  #   R's cairo DLL links against XQuartz (/opt/X11/lib/libXrender), so on a Mac
  #   without XQuartz `capabilities("cairo")` still reports TRUE while
  #   cairo_pdf() fails to load and R falls back to the base pdf() device. That
  #   device has no glyph for U+2212 or U+2013 and substitutes a hyphen for
  #   both, silently, with only an mbcsToSbcs warning. In this figure set that
  #   corrupts the axis labels "FR_RD (biomass - density)" and the condition
  #   classes "VeryPoor-Poor". cairosvg carries its own cairo and does not
  #   depend on XQuartz, so converting the SVG we just wrote preserves them.
  pdf_route <- if (file.exists(pdf_path)) "cairo_pdf" else NA_character_
  if (!file.exists(pdf_path) && nzchar(CAIROSVG) && file.exists(svg_path)) {
    st <- suppressWarnings(system2(CAIROSVG, c(shQuote(svg_path), "-o", shQuote(pdf_path)),
                                   stdout = FALSE, stderr = FALSE))
    if (identical(st, 0L) && file.exists(pdf_path)) pdf_route <- "cairosvg"
  }

  # FALLBACK 2 -- the base device, which is lossy. Say so every time.
  if (!file.exists(pdf_path)) {
    ggplot2::ggsave(pdf_path, p, width = w, height = h)
    pdf_route <- "base pdf() -- LOSSY"
    say("  NOTE: cairo unavailable and cairosvg not found; %s.pdf used the default", name)
    say("        device, which substitutes hyphens for en-dashes and minus signs.")
    say("        The PNG and the SVG are correct; use one of those for review.")
  }
  say("  wrote %s (%.1f x %.1f in)  png/svg/pdf [pdf via %s]", name, w, h, pdf_route)
}

#' Console + file logging, so a run leaves an audit trail.
LOGFILE <- file.path(OUT_LOG, format(Sys.time(), "run_%Y%m%d_%H%M%S.log"))
say <- function(...) {
  msg <- paste0(format(Sys.time(), "[%H:%M:%S] "), sprintf(...))
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOGFILE, append = TRUE)
  invisible(NULL)
}

say("config loaded | data folder: %s", if (REPO_LAYOUT) "data/ beside R/ (repository layout)" else DATA)

# ---- Trait group and state names as Twardochleb et al. (2021) publish them, 25 September 2026 -------------
# One map for every figure and table that prints a trait state; the codes are the trait matrix's. TRAIT_GROUP_SHORT holds
# the shortenings the manuscript allows in figures and Results: 'flying strength', 'emergence synchrony',
# 'body size'. TRAIT_GROUP_ORDER is the order the figures list the groups (by name).
TRAIT_GROUP_LAB <- c(Flying = "adult flying strength", Disperse = "female dispersal", Synchro = "synchronization of emergence",
                     Emerge = "emergence season", Feed = "feeding style", Habit = "habit", Size = "maximum body size",
                     Respiration = "respiration mode", Rheophily = "rheophily", Thermal = "thermal preference", Voltinism = "voltinism")
TRAIT_GROUP_SHORT <- c(Flying = "flying strength", Disperse = "female dispersal", Synchro = "emergence synchrony",
                       Emerge = "emergence season", Feed = "feeding style", Habit = "habit", Size = "body size",
                       Respiration = "respiration mode", Rheophily = "rheophily", Thermal = "thermal preference", Voltinism = "voltinism")
TRAIT_GROUP_ORDER <- c("Flying", "Disperse", "Synchro", "Emerge", "Feed", "Habit", "Size", "Respiration", "Rheophily", "Thermal", "Voltinism")
TRAIT_STATE_LAB <- c(Flying_Strong = "strong", Flying_Weak = "weak", Disperse_High = "high", Disperse_Low = "low",
                     Synchro_Well = "well synchronized", Synchro_Poorly = "poorly synchronized",
                     Emerge_Winter = "winter", Emerge_Spring = "spring", Emerge_Summer = "summer", Emerge_Fall = "fall",
                     Feed_SH = "shredder", Feed_CF = "collector-filterer", Feed_CG = "collector-gatherer", Feed_PR = "predator",
                     Feed_HB = "herbivore", Feed_PA = "parasite",
                     Habit_Burrower = "burrower", Habit_Climber = "climber", Habit_Clinger = "clinger", Habit_Crawler = "crawler",
                     Habit_Planktonic = "planktonic", Habit_Skater = "skater", Habit_Sprawler = "sprawler", Habit_Swimmer = "swimmer",
                     Size_Small = "small", Size_Medium = "medium", Size_Large = "large",
                     Respiration_Gills = "gills", Respiration_Tegument = "tegument", Respiration_Plastron = "plastron or spiracle",
                     Rheophily_Depositional = "depositional", Rheophily_Erosional = "erosional", Rheophily_Both = "both",
                     Thermal_ColdSteno = "cold stenothermal", Thermal_ColdEury = "cold-cool eurythermal",
                     Thermal_CoolEury = "cool-warm eurythermal", Thermal_WarmEury = "warm eurythermal", Thermal_HotEury = "hot eurythermal",
                     Voltinism_Semivoltine = "semivoltine", Voltinism_Univoltine = "univoltine", Voltinism_Multivoltine = "multivoltine")
stopifnot(length(TRAIT_STATE_LAB) == 41, setequal(names(TRAIT_GROUP_LAB), TRAIT_GROUP_ORDER))
