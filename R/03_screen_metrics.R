# =============================================================================
# 03_screen_metrics.R — within-region Spearman screen, BH within family
# =============================================================================
# WHY THIS STEP EXISTS
#   This is the analysis the whole paper rests on: each metric against B-IBI,
#   within each region, Spearman rank correlation, Benjamini-Hochberg correction
#   applied within each metric family within each region (Methods 2.5; see
#   add_bh_within_family in 00_config.R for why the family is region-specific).
#
#   v2 adds the concentration family from step 02 as a FIFTH family, so no
#   published FDR value moves.
#
# SCREENING RULE (Methods 2.5, unchanged)
#   strong       : |rho| >= 0.30 AND FDR-adjusted p <= 0.10
#   exploratory  : |rho| >= 0.30 AND raw p <= 0.05 but FDR > 0.10
#
# OUTPUT (step 07 builds Table 2 (legacy Tables 3a, 3b) and Table S4 part A (legacy Table S10) from these)
#   screen : one row per metric x region
#   output/tables/03_screen_all_metrics.csv
#   output/tables/03_screen_responsive_only.csv
#   output/tables/03_bh_scope_sensitivity.csv
# =============================================================================

say("STEP 03 -- screen metrics against B-IBI")

families <- c(METRIC_FAMILIES, list(concentration = CONCENTRATION_FAMILY))

dat_screen <- dat %>% left_join(conc, by = "SiteID")

screen <- map_dfr(names(families), function(fam) {
  map_dfr(families[[fam]], function(m) {
    if (!m %in% names(dat_screen)) {
      warning("metric not found, skipped: ", m); return(NULL)
    }
    map_dfr(REGIONS, function(rg) {
      s <- dat_screen[dat_screen$Region == rg, ]
      spearman_one(s$BIBI, s[[m]], rg, m) %>% mutate(Family = fam)
    })
  })
}) %>%
  add_bh_within_family(scope = BH_SCOPE) %>%
  add_tier() %>%
  mutate(
    # 27 Sep 2026: the weighting is the metric dictionary's, not the column suffix's, so presence-based
    # metrics read "invariant", the divergences "difference", and the tables agree with Table S2 part C (legacy Table S7).
    Weighting = vapply(Metric, function(m) mdict(m, "currency"), character(1)),
    Region = factor(Region, levels = REGIONS)
  ) %>%
  arrange(Family, Metric, Region) %>%
  select(Family, Metric, Weighting, Region, n, rho, p_raw, p_FDR, Tier)

write_csv(screen, file.path(OUT_TAB, "03_screen_all_metrics.csv"))
write_csv(screen %>% filter(Tier %in% c("strong", "exploratory")),
          file.path(OUT_TAB, "03_screen_responsive_only.csv"))

# ---- Guard: adding the new family must not have moved anything ---------------
# Re-check a few published FDR values. If these drift, the family grouping broke.
fdr_check <- tribble(
  ~Metric,      ~Region,          ~fdr_published,
  "FR_RD_bio",  "Piedmont",        0.0751,
  "FR_RD_den",  "Coastal Plain",   0.0020,
  "FDiv_den",   "Piedmont",        0.0497,
  "Shannon",    "Piedmont",        0.0006,
  "Total_Biomass","Piedmont",      0.0283
) %>%
  left_join(screen %>% select(Metric, Region, p_FDR) %>%
              mutate(Region = as.character(Region)), by = c("Metric","Region")) %>%
  mutate(ok = abs(p_FDR - fdr_published) < 0.002)

if (BH_SCOPE == "within_region" && !nzchar(Sys.getenv("CH2_BIBI_MEDIAN"))) {
  if (!all(fdr_check$ok)) {
    print(as.data.frame(fdr_check))
    stop("STEP 03 FAILED: published FDR values moved under v1 BH scope.\n",
         "  Either the concentration family leaked into an existing family, or an\n",
         "  input is wrong. Do not proceed.")
  }
  say("  provenance check passed: the new concentration family did not move any published FDR")
} else if (nzchar(Sys.getenv("CH2_BIBI_MEDIAN"))) {
  say("  reference FDR values were computed on the observed index; on the rarefied-median basis they are not enforced")
} else {
  say("  BH_SCOPE = '%s': reference FDR values intentionally not reproduced (see 00_config.R)", BH_SCOPE)
}

n_strong <- screen %>% filter(Tier == "strong") %>% count(Region)
say("  strong metrics by region (BH_SCOPE = %s): %s", BH_SCOPE,
    paste(sprintf("%s=%d", n_strong$Region, n_strong$n), collapse = "  "))

# ---- The other BH scope, computed and reported but NOT adopted ---------------
# Methods 2.5 literally says the correction is applied "within each
# metric-hypothesis family ... across regions". v1 applied it within family AND
# within region. This block shows exactly what would change if the wording were
# taken at face value; the paper keeps BH_SCOPE = "within_region".
alt_scope <- setdiff(c("within_region","across_regions"), BH_SCOPE)
screen_alt <- screen %>% select(Family, Metric, Weighting, Region, n, rho, p_raw) %>%
  add_bh_within_family(scope = alt_scope) %>% add_tier()

tier_delta <- screen %>% select(Metric, Region, Tier_current = Tier, FDR_current = p_FDR) %>%
  left_join(screen_alt %>% select(Metric, Region, Tier_alt = Tier, FDR_alt = p_FDR),
            by = c("Metric","Region")) %>%
  filter(Tier_current != Tier_alt) %>%
  mutate(BH_scope_current = BH_SCOPE, BH_scope_alt = alt_scope)

write_csv(tier_delta, file.path(OUT_TAB, "03_bh_scope_sensitivity.csv"))

if (nrow(tier_delta)) {
  say("  BH scope sensitivity: %d metric x region results change tier under '%s'",
      nrow(tier_delta), alt_scope)
  for (i in seq_len(nrow(tier_delta)))
    say("    %-18s %-15s %s -> %s   (FDR %.4f -> %.4f)",
        tier_delta$Metric[i], tier_delta$Region[i], tier_delta$Tier_current[i],
        tier_delta$Tier_alt[i], tier_delta$FDR_current[i], tier_delta$FDR_alt[i])
  say("  This is a decision, not a bug. BH_SCOPE is set in 00_config.R; the alternative is reported for comparison only.")
} else {
  say("  BH scope sensitivity: no result changes tier under '%s'", alt_scope)
}
say("STEP 03 complete")
