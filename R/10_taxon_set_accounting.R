# =============================================================================
# 10_taxon_set_accounting.R — define the analyzed set, account for the rest
# =============================================================================
# WHY THIS STEP EXISTS
#   The analysis rests on the trait-eligible matrix (the taxa with a trait
#   record). The community table holds more taxa than that. The excluded taxa
#   are not noise: they are the non-insect macroinvertebrates for which the
#   continent-scale trait database has no entry -- isopods, amphipods,
#   gastropods, bivalves, mites. They carry real density and real biomass. The
#   counts are computed below, never typed.
#
#   This step:
#     1. States the analyzed set and the excluded set explicitly, with numbers.
#     2. Reports trait coverage per site against the standard convention.
#     3. Tests whether each responsive association survives the convention's threshold.
#     4. Records each family's share of density and biomass (section 5).
#
# THE CONVENTION
#   Pakeman & Quested (2007, Applied Vegetation Science 10:91-96) established
#   what is now the standard criterion for trait-based community analysis: the
#   taxa with trait data should account for at least 80% of community abundance
#   in each sample. Below that, community-weighted indices become sensitive to
#   which taxa are missing rather than to the gradient. Reporting coverage and
#   testing against the 80% threshold is what the trait literature expects.
#
# WHY IT MATTERS HERE, SPECIFICALLY
#   (a) Trait coverage is the share of community density or biomass carried by
#       the trait-eligible taxa (the columns of den_mat and bio_mat), computed
#       below for every site, never typed.
#
#   (b) The excluded fraction differs among regions, and the paper's argument
#       is a regional comparison, so the exclusion is stated by region.
#
# TOTALS UNDER TWO TAXON SETS
#   Total_Density and Total_Biomass in the site-summary file are row sums of the
#   full community matrix (every taxon). Every trait-based metric, and every
#   multivariate test, is computed on the trait-eligible matrix. Section 3 of
#   this step writes the correlation of each total with condition under both
#   taxon sets, so the difference is visible.
#
# OUTPUT
#   output/tables/10_site_trait_coverage.csv       repository table trait_coverage_by_site.csv (legacy Table S4); also read by step 18
#   output/tables/10_excluded_taxa_accounting.csv  log record; no later step reads it
#   output/tables/10_totals_two_definitions.csv    log record; no later step reads it
#   output/tables/10_coverage80_sensitivity.csv    supplement Table S3 part A (legacy Table S5), joined with step 12 by step 18
#   output/tables/10_family_density_share.csv      repository table family_density_biomass_shares.csv (Discussion 4.5)
# =============================================================================

say("STEP 10 -- taxon-set accounting and trait coverage")

COVERAGE_MIN <- 0.80   # Pakeman & Quested (2007)

raw <- read_csv(file.path(DATA, "SiteTaxa_Density_Biomass_Transformed.csv"),
                show_col_types = FALSE) %>%
  mutate(SiteID = std_site(SiteID), Region = as_region(Typology),
         den = coalesce(Median_Density_per_m2, 0),
         bio = coalesce(Median_Biomass_per_m2, 0),
         analysed = TAXANAME %in% colnames(den_mat))

# ---- 1. Per-site coverage ----------------------------------------------------
cov <- raw %>% group_by(SiteID, Region) %>%
  summarise(n_taxa_total = n_distinct(TAXANAME),
            n_taxa_analyzed = n_distinct(TAXANAME[analysed]),
            den_coverage = sum(den[analysed]) / sum(den),
            bio_coverage = sum(bio[analysed]) / sum(bio), .groups = "drop") %>%
  mutate(meets_80_density = den_coverage >= COVERAGE_MIN,
         meets_80_biomass = bio_coverage >= COVERAGE_MIN,
         meets_80_both    = meets_80_density & meets_80_biomass)
write_csv(cov, file.path(OUT_TAB, "10_site_trait_coverage.csv"))

cov_sum <- cov %>% group_by(Region) %>%
  summarise(den_median = median(den_coverage), den_min = min(den_coverage),
            bio_median = median(bio_coverage), bio_min = min(bio_coverage),
            n_below_80 = sum(!meets_80_both), n = n(), .groups = "drop")
say("  trait coverage by region (share of community carried by the %d analyzed taxa):", ncol(den_mat))
for (i in seq_len(nrow(cov_sum))) {
  r <- cov_sum[i, ]
  say("    %-15s density %.1f%% median (min %.1f%%) | biomass %.1f%% median (min %.1f%%) | %d of %d sites below 80%%",
      r$Region, 100*r$den_median, 100*r$den_min, 100*r$bio_median, 100*r$bio_min, r$n_below_80, r$n)
}
say("  coverage is the share of site density or biomass carried by the trait-eligible taxa")

# ---- 2. The excluded component ----------------------------------------------
# Broad functional identity of the excluded groups, assigned from taxonomy.
# These are stated so the gap can be named in the Discussion rather than left
# as an unexplained exclusion. None of them is a crayfish.
group_of <- function(x) case_when(
  grepl("Lirceus|Asellus|Caecidotea", x)                        ~ "Isopoda (detritivore/shredder)",
  grepl("Crangonyx|Gammarus|Hyalella|Amphipoda", x)             ~ "Amphipoda (detritivore/shredder)",
  grepl("Physidae|Lymnaeidae|Pleuroceridae|Viviparidae|Ancylidae|Bithynia|Gastropoda|Planorbidae", x)
                                                                ~ "Gastropoda (scraper/grazer)",
  grepl("Corbicula|Sphaeriidae|Bivalvia|Pisidium", x)           ~ "Bivalvia (filter-feeder)",
  grepl("Acari|Hydracarina", x)                                 ~ "Acari (predator/parasite)",
  grepl("Collembola|Chilopoda|Aphididae|Curculionidae|Chrysomelidae", x)
                                                                ~ "Terrestrial/incidental",
  TRUE                                                          ~ "Other (no trait record)"
)

excluded <- raw %>% filter(!analysed) %>%
  mutate(FunctionalGroup = group_of(TAXANAME)) %>%
  group_by(Region, FunctionalGroup) %>%
  summarise(n_taxa = n_distinct(TAXANAME), n_sites = n_distinct(SiteID),
            den = sum(den), bio = sum(bio), .groups = "drop") %>%
  left_join(raw %>% group_by(Region) %>% summarise(den_tot = sum(den), bio_tot = sum(bio), .groups = "drop"),
            by = "Region") %>%
  mutate(pct_regional_density = 100 * den / den_tot,
         pct_regional_biomass = 100 * bio / bio_tot) %>%
  select(Region, FunctionalGroup, n_taxa, n_sites,
         pct_regional_density, pct_regional_biomass) %>%
  arrange(Region, desc(pct_regional_biomass))
write_csv(excluded, file.path(OUT_TAB, "10_excluded_taxa_accounting.csv"))

reg_excl <- raw %>% group_by(Region) %>%
  summarise(pct_den = 100*sum(den[!analysed])/sum(den),
            pct_bio = 100*sum(bio[!analysed])/sum(bio),
            n_taxa = n_distinct(TAXANAME[!analysed]), .groups = "drop")
say("  excluded component, share of the regional total:")
for (i in seq_len(nrow(reg_excl))) {
  r <- reg_excl[i, ]
  say("    %-15s %.1f%% of density, %.1f%% of biomass, across %d taxa",
      r$Region, r$pct_den, r$pct_bio, r$n_taxa)
}
say("  the excluded share is stated by region because the paper is a regional comparison.")

# ---- 3. Totals under both definitions ----------------------------------------
two_def <- dat %>%
  mutate(TE_den = rowSums(den_mat[SiteID, , drop = FALSE]),
         TE_bio = rowSums(bio_mat[SiteID, , drop = FALSE])) %>%
  group_by(Region) %>%
  summarise(
    den_trait_set_rho = cor(BIBI, TE_den, method = "spearman"),
    den_all_taxa_rho = cor(BIBI, Total_Density, method = "spearman"),
    bio_trait_set_rho = cor(BIBI, TE_bio, method = "spearman"),
    bio_all_taxa_rho = cor(BIBI, Total_Biomass, method = "spearman"),
    .groups = "drop") %>%
  mutate(n_taxa_trait_set = ncol(den_mat), n_taxa_all = n_distinct(raw$TAXANAME))
write_csv(two_def, file.path(OUT_TAB, "10_totals_two_definitions.csv"))

say("  standing-stock totals under the two taxon sets (trait set, %d taxa; every taxon, %d):", ncol(den_mat), n_distinct(raw$TAXANAME))
for (i in seq_len(nrow(two_def))) {
  r <- two_def[i, ]
  say("    %-15s density  trait set %+.3f   all taxa %+.3f", r$Region, r$den_trait_set_rho, r$den_all_taxa_rho)
  say("                    biomass  trait set %+.3f   all taxa %+.3f", r$bio_trait_set_rho, r$bio_all_taxa_rho)
}

# ---- 4. Does each responsive association survive the 80% rule? ----------------
keep_sites <- cov %>% filter(meets_80_both) %>% pull(SiteID)
dat80 <- dat %>% left_join(conc, by = "SiteID") %>% filter(SiteID %in% keep_sites)

# 28 Sep 2026: the associations re-tested are the screen's responsive ones (strong or exploratory tier in
# step 03's 03_screen_all_metrics.csv), metric by region, as step 12 does. A typed list of 11 "headline" metrics had missed
# associations that meet the criterion (richness, HBI by biomass, Gini by biomass, N50, FRic, FR_SF) and re-tested some that do not.
resp80 <- read_csv(file.path(OUT_TAB, "03_screen_all_metrics.csv"), show_col_types = FALSE) %>%
  filter(Tier %in% c("strong", "exploratory")) %>% select(Metric, Region)
stopifnot(nrow(resp80) > 0)
sens80 <- purrr::pmap_dfr(resp80, function(Metric, Region) {
  m <- Metric; rg <- Region
  a <- dat %>% left_join(conc, by = "SiteID") %>% filter(Region == rg)
  b <- dat80 %>% filter(Region == rg)
  tibble(Metric = m, Region = rg,
         n_all = nrow(a), rho_all = cor(a$BIBI, a[[m]], method = "spearman", use = "complete.obs"),
         n_80  = nrow(b), rho_80  = if (nrow(b) >= 8) cor(b$BIBI, b[[m]], method = "spearman", use = "complete.obs") else NA_real_)
}) %>% mutate(shift = rho_80 - rho_all,
               sign_change = sign(rho_80) != sign(rho_all))
write_csv(sens80, file.path(OUT_TAB, "10_coverage80_sensitivity.csv"))

say("  80%%-rule sensitivity: %d of %d sites retained", length(keep_sites), nrow(cov))
big <- sens80 %>% filter(!is.na(shift), abs(shift) > 0.15 | sign_change)
if (nrow(big)) {
  say("  metrics shifting by more than 0.15 when low-coverage sites are dropped:")
  for (i in seq_len(nrow(big)))
    say("    %-14s %-15s all sites %+.3f (n=%d) -> >=80%% coverage %+.3f (n=%d)",
        big$Metric[i], big$Region[i], big$rho_all[i], big$n_all[i], big$rho_80[i], big$n_80[i])
} else {
  say("  every responsive association is stable under the 80%% criterion -- no shift exceeds 0.15.")
}
# ---- 5. Family shares of density and biomass (8 Oct 2026; Discussion 4.5) ------------------
# The Discussion says Chironomidae made up the largest share of density in every region. This table is the record behind
# that sentence: each family's share of the regional density and biomass (sum over streams of the site median), and the
# number of streams in which the family has the largest share. Every taxon counts, not only the trait set.
fam <- raw %>% mutate(Family = toupper(trimws(Family)))
fam_share <- fam %>% group_by(Region, Family) %>% summarise(den = sum(den), bio = sum(bio), .groups = "drop") %>%
  group_by(Region) %>% mutate(share_density = den / sum(den), share_biomass = bio / sum(bio)) %>% ungroup()
fam_top <- fam %>% group_by(Region, SiteID, Family) %>% summarise(den = sum(den), .groups = "drop") %>%
  group_by(Region, SiteID) %>% slice_max(den, n = 1, with_ties = FALSE) %>% ungroup() %>%
  count(Region, Family, name = "streams_largest_density_share")
fam_share <- fam_share %>% left_join(fam_top, by = c("Region", "Family")) %>%
  mutate(streams_largest_density_share = coalesce(streams_largest_density_share, 0L),
         rank_density = ave(-share_density, Region, FUN = function(v) rank(v, ties.method = "min"))) %>%
  select(Region, Family, rank_density, share_density, share_biomass, streams_largest_density_share) %>%
  arrange(factor(Region, levels = REGIONS), rank_density)
write_csv(fam_share, file.path(OUT_TAB, "10_family_density_share.csv"))
for (r in REGIONS) { f1 <- fam_share %>% filter(Region == r, rank_density == 1)
  say("  largest family share of density, %s: %s (%.1f%% of density; largest in %d of %d streams)", r, f1$Family,
      100 * f1$share_density, f1$streams_largest_density_share, n_distinct(fam$SiteID[fam$Region == r])) }

say("STEP 10 complete")
