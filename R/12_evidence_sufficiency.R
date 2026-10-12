# =============================================================================
# 12_evidence_sufficiency.R — is a responsive association a trait signal or a coverage artifact?
# =============================================================================
# WHY THIS STEP EXISTS
#   A taxon with no trait record contributes nothing to a trait-based metric, so a
#   metric can track condition because the SHARE of the community with trait data
#   tracks condition. That share is not random with respect to condition: sites in
#   better condition hold proportionally more insects and fewer of the non-insect
#   taxa the trait database has no record for.
#
# THE TEST
#   For every screen association at the strong or exploratory tier: the share of
#   the site's density (or biomass, for a biomass-weighted metric) in taxa with a
#   trait record; the metric's correlation with that share; and the metric's
#   partial rank correlation with condition after the share is held constant. A
#   trait signal survives; a coverage artifact does not. The correlation of
#   coverage itself with condition, in each region, is written first.
#
# OUTPUT
#   output/tables/12_coverage_vs_condition_by_region.csv   log record; no later step reads it
#   output/tables/12_coverage_confound_test.csv    joined with step 10's 80% subset by step 18 into supplement
#                                                  Table S3 part A (legacy Table S5; formerly its own Table S21)
# =============================================================================

say("STEP 12 -- evidence sufficiency: trait coverage partialled out of every responsive association")
# Until 24 September 2026 this step tested legacy community-weighted means from the retired binary file, none of which the
# paper reports. It now runs the same partial rank correlation on the associations the paper does
# report: every screen result at the strong or exploratory tier (03_screen_all_metrics.csv), with the share of site
# density (or biomass, for a biomass-weighted metric) in taxa that have a trait record held constant (the coverage is
# computed here from the same community table step 10 reads, so this step does not depend on step 10).
raw <- read_csv(file.path(DATA, "SiteTaxa_Density_Biomass_Transformed.csv"),
                show_col_types = FALSE) %>%
  mutate(SiteID = std_site(SiteID), den = coalesce(Median_Density_per_m2, 0),
         bio = coalesce(Median_Biomass_per_m2, 0), an = TAXANAME %in% colnames(den_mat))
cov <- raw %>% group_by(SiteID) %>%
  summarise(den_cov = sum(den[an])/sum(den), bio_cov = sum(bio[an])/sum(bio), .groups = "drop")
d <- dat %>% left_join(conc, by = "SiteID") %>% inner_join(cov, by = "SiteID")
say("  trait coverage vs condition, by region (the reason the confound matters):")
cov_reg <- purrr::map_dfr(REGIONS, function(rg) {
  s <- d %>% filter(Region == rg)
  tibble(Region = rg, rho_density_coverage_vs_index = cor(s$BIBI, s$den_cov, method = "spearman"),
         rho_biomass_coverage_vs_index = cor(s$BIBI, s$bio_cov, method = "spearman"),
         excluded_share_density_median = 100 * (1 - median(s$den_cov)), excluded_share_biomass_median = 100 * (1 - median(s$bio_cov)))
})
write_csv(cov_reg, file.path(OUT_TAB, "12_coverage_vs_condition_by_region.csv"))
for (i in seq_len(nrow(cov_reg))) say("    %-15s density coverage rho = %+.3f | biomass coverage rho = %+.3f", cov_reg$Region[i], cov_reg$rho_density_coverage_vs_index[i], cov_reg$rho_biomass_coverage_vs_index[i])
resp <- screen %>% filter(Tier %in% c("strong", "exploratory"))
res <- purrr::pmap_dfr(resp %>% select(Family, Metric, Weighting, Region, Tier, rho, p_FDR), function(Family, Metric, Weighting, Region, Tier, rho, p_FDR) {
  s  <- d %>% filter(Region == !!Region)
  cv <- if (Weighting == "Biomass") s$bio_cov else s$den_cov          # by the screen's weighting column, not the name suffix (27 Sep 2026)
  r_marg <- cor(s$BIBI, s[[Metric]], method = "spearman", use = "complete.obs")
  r_part <- cor(resid(lm(rank(s[[Metric]]) ~ rank(cv))), resid(lm(rank(s$BIBI) ~ rank(cv))))
  tibble(Family, Metric, Region, Tier, coverage_used = ifelse(Weighting == "Biomass", "biomass", "density"),
         r_metric_with_coverage = cor(s[[Metric]], cv, use = "complete.obs", method = "spearman"),
         rho_marginal = r_marg, rho_partial = r_part, pct_lost = 100 * (1 - abs(r_part) / abs(r_marg)),
         keeps_criterion_coefficient = abs(r_part) >= RHO_MIN,
         verdict = ifelse(abs(r_part) >= RHO_MIN, "keeps |rho| >= 0.30 with coverage held constant", "falls below |rho| 0.30 with coverage held constant"))
})
write_csv(res, file.path(OUT_TAB, "12_coverage_confound_test.csv"))
say("  %d responsive associations: %d keep |rho| >= %.2f with coverage held constant, %d fall below it", nrow(res), sum(res$keeps_criterion_coefficient), RHO_MIN, sum(!res$keeps_criterion_coefficient))
for (i in which(!res$keeps_criterion_coefficient)) say("    %-18s %-15s %-12s rho %+.3f -> %+.3f partial", res$Metric[i], res$Region[i], res$Tier[i], res$rho_marginal[i], res$rho_partial[i])
say("STEP 12 complete")
