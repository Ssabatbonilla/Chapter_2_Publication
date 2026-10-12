# =============================================================================
# 07_tables.R — screen and multivariate tables that the manuscript tables are built from
# =============================================================================
# WHY THIS STEP EXISTS
#   Turns the screen (step 03) and the multivariate tests (step 04) into the
#   tables the paper carries, with the field's names for every metric and
#   family. Step 18 numbers and places them:
#     07_T2_responsive_metrics.csv    every association in the strong or exploratory tier -> Table 2
#                                     (step 18 splits it into legacy Tables 3a and 3b; step 18b merges them)
#     07_S1_full_screen.csv           the full screen, all metrics and regions -> Table S4 part A (legacy S10)
#     07_S2_regional_distributions.csv  regional metric distributions       -> repository table
#                                     regional_metric_distributions.csv (legacy Table S13; not in the supplement)
#     07_S4_multivariate_tests.csv    multivariate tests, cross-region and within-region -> Table S5 (legacy S9)
#   The study-system table (Table S1 part A) is built by step 18. No land-use table is
#   built: no Chapter 2 analysis uses land use.
#
# OUTPUT: output/tables/07_*.csv
# =============================================================================

say("STEP 07 -- manuscript tables")

pretty <- c(
  Richness_alpha = "Taxonomic richness", Shannon = "Shannon diversity (H')",
  HBI_rel_abund  = "HBI (relative abundance)",
  HBI_bio        = "HBI (biomass)",
  Total_Density  = "Total density", Total_Biomass = "Total biomass",
  FRic_den = "Functional richness (FRic)", FEve_den = "Functional evenness (FEve)",
  FEve_bio = "Functional evenness (FEve)", FDiv_den = "Functional divergence (FDiv)",
  FDiv_bio = "Functional divergence (FDiv)", FDis_den = "Functional dispersion (FDis)",
  FDis_bio = "Functional dispersion (FDis)", RaoQ_den = "Rao's quadratic entropy",
  RaoQ_bio = "Rao's quadratic entropy", FR_RD_den = "Overlap redundancy (FR_RD)",
  FR_RD_bio = "Overlap redundancy (FR_RD)", FR_SF_den_log10 = "Space-filling redundancy (FR_SF)",
  FR_DB_RD = "Redundancy divergence", Gini_den = "Gini coefficient",
  Gini_bio = "Gini coefficient", N50_den = "Taxa making up half the total",
  N50_bio = "Taxa making up half the total", Gini_DB_DOM = "Concentration divergence"
)
property <- c(alpha = "Structure and tolerance", structure = "Totals",   # 25 Sep 2026: the keys stay (they are the BH families); the display names are the field's
              fd = "Functional diversity", redundancy = "Redundancy",
              concentration = "Inequality")

# ---- T2: the merged responsive-metrics table --------------------------------
# HBI weighted by areal density is the conventional relative-abundance HBI:
# density is a constant multiple of relative abundance within a site, so the two
# give identical rank correlations to 16 significant figures in every region.
# Reporting both inflated the structural family from four metrics to five. The
# relative-abundance form is the one Hilsenhoff (1988) defines, so it is kept and
# the density label is dropped. HBI_den is removed in 00_config.R and does not
# exist anywhere in the pipeline, so no filter is needed here.
t2 <- screen %>%
  filter(Tier %in% c("strong", "exploratory")) %>%
  mutate(Property = property[Family],
         Metric_label = ifelse(Metric %in% names(pretty), pretty[Metric], Metric),
         `Spearman rho` = round(rho, 3),
         `raw p` = signif(p_raw, 3),
         `FDR p` = signif(p_FDR, 3)) %>%
  select(Property, Metric = Metric_label, Weighting, Region, n,
         `Spearman rho`, `raw p`, `FDR p`, Tier) %>%
  arrange(factor(Property, levels = property), Metric, Region)
write_csv(t2, file.path(OUT_TAB, "07_T2_responsive_metrics.csv"))
say("  T2 responsive metrics: %d rows (the source of Table 2, via steps 18 and 18b)", nrow(t2))

# ---- S1: complete screen ----------------------------------------------------
# NOTE ON SCOPE. An earlier version of this table also carried 82 community-weighted
# trait means, which belong to a companion paper. This table deliberately contains
# only the screened metrics of METRIC_FAMILIES (19) plus the 5 concentration
# metrics. If you want the trait-mean screen published here, add
# it back knowingly -- do not let it arrive by default.
s1 <- screen %>%
  mutate(Property = property[Family],
         Metric_label = ifelse(Metric %in% names(pretty), pretty[Metric], Metric)) %>%
  select(Property, Family, Metric, Metric_label, Weighting, Region, n,
         rho, p_raw, p_FDR, Tier) %>%
  arrange(Family, Metric, Region)
write_csv(s1, file.path(OUT_TAB, "07_S1_full_screen.csv"))
say("  S1 full screen: %d metrics x %d regions = %d rows (trait-mean screen deliberately excluded)",
    n_distinct(s1$Metric), length(REGIONS), nrow(s1))

# ---- S2: regional distributions ---------------------------------------------
s2 <- dat %>% left_join(conc, by = "SiteID") %>%
  select(Region, any_of(c(unlist(METRIC_FAMILIES, use.names = FALSE), CONCENTRATION_FAMILY))) %>%
  pivot_longer(-Region, names_to = "Metric", values_to = "value") %>%
  filter(is.finite(value)) %>%
  group_by(Metric, Region) %>%
  summarise(n = n(), median = median(value), min = min(value), max = max(value),
            q25 = quantile(value, .25), q75 = quantile(value, .75), .groups = "drop") %>%
  mutate(`median [min-max]` = sprintf("%.3f [%.3f-%.3f]", median, min, max)) %>%
  arrange(Metric, Region)
write_csv(s2, file.path(OUT_TAB, "07_S2_regional_distributions.csv"))
say("  S2 regional distributions: %d rows", nrow(s2))

# S3 (land-use covariate models, from a former step 05) is no longer built: no Chapter 2 analysis uses land use
# (24 Sep 2026). Step 05 is not part of this repository.

# S4 combines both of step 04's models. Section 3.1 reports the cross-region
# model in its first paragraph and the within-region models in its second, so
# one table carrying both is what the section cites, and it makes visible
# that the same transformation was applied to both.
s4 <- bind_rows(
  read_csv(file.path(OUT_TAB, "04_cross_region_multivariate.csv"), show_col_types = FALSE) %>%
    mutate(Scope = "cross-region", .before = 1),
  read_csv(file.path(OUT_TAB, "04_within_region_multivariate.csv"), show_col_types = FALSE) %>%
    mutate(Scope = "within-region", .before = 1)
) %>% arrange(Scope, Transform, Weighting, Region, Test, Term)
write_csv(s4, file.path(OUT_TAB, "07_S4_multivariate_tests.csv"))
say("  S4 combined from step 04 (%d rows: %d cross-region, %d within-region)",
    nrow(s4), sum(s4$Scope == "cross-region"), sum(s4$Scope == "within-region"))
say("STEP 07 complete")
