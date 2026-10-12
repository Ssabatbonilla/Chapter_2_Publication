# =============================================================================
# 11_taxon_convention_and_traitspace_quality.R
# =============================================================================
# PART 1 -- THE TWO-TIER TAXON CONVENTION (decided 2026-08-26)
#
#   Community-level analyses -- composition, and the standing-stock totals --
#   use EVERY taxon that was quantitatively sampled. Trait-level analyses use
#   the subset with resolved trait profiles. The two tiers are different by
#   necessity, and the paper states which is which rather than blurring them.
#
#     COMMUNITY SET  = all taxa collected, minus those not sampled
#                      quantitatively (Decapoda, already removed upstream) and
#                      minus terrestrial incidentals that entered by drift.
#     TRAIT SET      = the taxa with a resolved trait profile (the columns of
#                      den_mat; the count is written to 11_community_vs_trait_set.csv).
#
#   Why this is the right split. A total is a statement about the community, so
#   it should be computed on the community. A trait-space property is a
#   statement about the trait-resolved fraction, and cannot honestly be
#   computed on anything else. Forcing both onto one taxon set means either
#   discarding real biomass from the totals or inventing trait values.
#
#   The cost is that the standing-stock totals and the trait metrics describe different assemblages. That cost is
#   paid openly, in Methods and in the Discussion, and it converts into the
#   paper's forward-looking claim: extending trait coverage to non-insect
#   macroinvertebrates is a prerequisite for a complete functional assessment.
#
#   TERRESTRIAL EXCLUSIONS are deliberately minimal (their share of biomass is logged below). The
#   aquatic non-insects (Oligochaeta, water mites, Turbellaria) are real benthic
#   fauna, quantitatively sampled, and STAY in the community set even though
#   they have no trait record. They are exactly the gap the Discussion names.
#
# PART 2 -- IS THE TRAIT SPACE ADEQUATELY REPRESENTED?
#
#   A trait-space ordination plots two PCoA axes, which carry a small share of
#   the positive eigenvalue variance (11_traitspace_dimensionality.csv). That
#   number alone does not say whether the picture is
#   faithful, because what matters is whether distances in the reduced space
#   reproduce the trait dissimilarities the metrics are computed from.
#
#   Maire et al. (2015, Global Ecol. Biogeogr. 24:728-740) give the standard
#   test: the mean squared deviation (mSD) between the original trait distance
#   and the Euclidean distance in k dimensions, both scaled to a maximum of 1.
#   This step runs it, with the mean absolute deviation and the correlation
#   between the two distances beside it.
#
#   THE ANSWER MATTERS FOR HOW THE FIGURE IS FRAMED, NOT FOR THE METRICS.
#   FD::dbFD computes RaoQ and FDis from the Gower distance and FDiv, FEve and
#   FRic from up to ten corrected principal coordinates (step 13), not from the
#   two plotted axes. So the coefficients are unaffected. What is affected is
#   how much inferential weight an ordination picture can carry.
#
# OUTPUT
#   output/tables/11_community_vs_trait_set.csv          read by step 18
#   output/tables/11_community_set_standing_stock.csv    log record; no later step reads it
#   output/tables/11_traitspace_dimensionality.csv       repository table traitspace_dimensionality.csv (legacy Table S6)
#   output/tables/11_traitspace_taxon_coordinates.csv, 11_traitspace_axes.csv, 11_traitspace_site_scores.csv
#                                                        read by step 06 for fig4_traitspace_binned (drawn, not in the paper)
#   output/tables/11_trait_group_correlations.csv        log record
#   output/tables/11_parent_rows_resolved.csv            log record
#   output/tables/11_fric_vs_richness.csv                supplement Table S4 part C (legacy Table S25)
# =============================================================================

say("STEP 11 -- taxon convention and trait-space quality")

# Terrestrial / drift taxa: collected, but not part of the benthic assemblage
# the design targets. Keep this list short and defensible.
TERRESTRIAL <- c("Chilopoda","Diplopoda","Aphididae","Curculionidae",
                 "Chrysomelidae","Staphylinidae","Thysanoptera","Formicidae")

raw <- read_csv(file.path(DATA, "SiteTaxa_Density_Biomass_Transformed.csv"),
                show_col_types = FALSE) %>%
  mutate(SiteID = std_site(SiteID), Region = as_region(Typology),
         den = coalesce(Median_Density_per_m2, 0),
         bio = coalesce(Median_Biomass_per_m2, 0),
         in_trait_set = TAXANAME %in% colnames(den_mat),
         terrestrial  = TAXANAME %in% TERRESTRIAL,
         in_community_set = !terrestrial)

# ---- 1. Define and describe the two sets ------------------------------------
sets <- raw %>%
  summarise(
    n_collected        = n_distinct(TAXANAME),
    n_terrestrial      = n_distinct(TAXANAME[terrestrial]),
    n_community_set    = n_distinct(TAXANAME[in_community_set]),
    n_trait_set        = n_distinct(TAXANAME[in_trait_set]),
    n_community_no_traits = n_distinct(TAXANAME[in_community_set & !in_trait_set]),
    pct_bio_terrestrial   = 100*sum(bio[terrestrial])/sum(bio),
    pct_bio_no_traits     = 100*sum(bio[in_community_set & !in_trait_set])/sum(bio),
    pct_den_no_traits     = 100*sum(den[in_community_set & !in_trait_set])/sum(den)
  )
write_csv(sets, file.path(OUT_TAB, "11_community_vs_trait_set.csv"))
say("  collected %d taxa | terrestrial excluded %d (%.2f%% of biomass) | COMMUNITY SET %d | TRAIT SET %d",
    sets$n_collected, sets$n_terrestrial, sets$pct_bio_terrestrial,
    sets$n_community_set, sets$n_trait_set)
say("  in the community set but without trait data: %d taxa, %.1f%% of density, %.1f%% of biomass",
    sets$n_community_no_traits, sets$pct_den_no_traits, sets$pct_bio_no_traits)

# ---- 2. Standing-stock totals on the COMMUNITY set --------------------------
tot <- raw %>% filter(in_community_set) %>%
  group_by(SiteID, Region) %>%
  summarise(Total_Density_community = sum(den),
            Total_Biomass_community = sum(bio), .groups = "drop") %>%
  left_join(dat %>% select(SiteID, BIBI), by = "SiteID")

ss <- map_dfr(c("Total_Density_community","Total_Biomass_community"), function(m)
  map_dfr(REGIONS, function(rg) {
    s <- tot %>% filter(Region == rg)
    ct <- suppressWarnings(cor.test(s$BIBI, s[[m]], method = "spearman", exact = FALSE))
    tibble(Metric = m, Region = rg, n = nrow(s),
           rho = unname(ct$estimate), p_raw = ct$p.value)
  })) %>%
  group_by(Region) %>% mutate(p_FDR = p.adjust(p_raw, method = "BH")) %>% ungroup() %>%
  mutate(Tier = ifelse(abs(rho) >= RHO_MIN & p_FDR <= FDR_MAX, "strong", "not responsive"))
write_csv(ss, file.path(OUT_TAB, "11_community_set_standing_stock.csv"))

say("  standing stock on the COMMUNITY set (every taxon except terrestrial incidentals):")
for (i in seq_len(nrow(ss)))
  say("    %-25s %-15s rho=%+.3f  FDR=%.4f  %s",
      ss$Metric[i], ss$Region[i], ss$rho[i], ss$p_FDR[i],
      ifelse(ss$Tier[i] == "strong", "STRONG", ""))

# ---- 3. Trait-space dimensionality (Maire et al. 2015) ----------------------
suppressPackageStartupMessages({library(FD); library(ade4)})
tr <- readRDS(IN$traits) %>% filter(TAXANAME %in% colnames(den_mat)) %>% as.data.frame()
rownames(tr) <- tr$TAXANAME; tr$TAXANAME <- NULL; tr$Genus <- NULL   # Twardochleb fuzzy file since 22 Sep 2026
tr <- tr[, colSums(tr, na.rm = TRUE) > 0, drop = FALSE]              # NA kept: an empty group is missing
tr <- as.data.frame(scale(tr))                    # dbFD uses stand.x = TRUE

g   <- suppressWarnings(gowdis(tr))
gm  <- as.matrix(g); ref <- gm[lower.tri(gm)]
pco <- ade4::dudi.pco(suppressWarnings(cailliez(g)), scannf = FALSE, nf = 30)
pos <- pco$eig[pco$eig > 0]; cum <- cumsum(pos)/sum(pos)

dims <- map_dfr(c(2,3,4,5,6,8,10,15,20), function(k) {
  co <- as.matrix(pco$li[, 1:k, drop = FALSE])
  de <- as.matrix(dist(co)); e_raw <- de[lower.tri(de)]
  e  <- e_raw * (mean(ref)/mean(e_raw))
  # Maire et al. (2015, p. 731): the mean squared deviation between the initial distance and the distance in the
  # k-axis space, both scaled to the same maximum (here 1), so mSD runs from 0 to 1 for a Gower distance
  mSD <- mean((ref / max(ref) - e_raw / max(e_raw))^2)
  tibble(n_axes = k,
         cum_variance_pct = 100*cum[k],
         mSD_maire = mSD,
         mAD = mean(abs(ref - e)),
         mAD_pct_of_mean_distance = 100*mean(abs(ref - e))/mean(ref),
         pearson_r = cor(ref, e),
         spearman_r = cor(ref, e, method = "spearman"))
})
write_csv(dims, file.path(OUT_TAB, "11_traitspace_dimensionality.csv"))

# ---- 3b. The trait-space ordination read by step 06 (27 Sep 2026) ----------------------------------------------------
# Step 06 draws fig4_traitspace_binned from these files; the figure is not in this paper (it was supplement Figure S1
# until 28 Sep 2026; the paper's Figure S1 is the chironomid-removal figure). The same principal coordinate analysis as
# above (the current fuzzy matrix, the taxa of the trait set, Cailliez correction) is written here: the taxon coordinates on the first ten axes, each axis's
# share of the positive eigenvalues, and each site's score as its abundance-weighted centroid of the taxa present,
# weighted by median density and by median biomass.
tx_co <- as.data.frame(pco$li[, 1:10, drop = FALSE]); names(tx_co) <- paste0("A", seq_len(ncol(tx_co)))
tx_co <- tibble(TAXANAME = rownames(pco$li)) %>% bind_cols(tx_co)
stopifnot(setequal(tx_co$TAXANAME, colnames(den_mat)))
write_csv(tx_co, file.path(OUT_TAB, "11_traitspace_taxon_coordinates.csv"))
write_csv(tibble(axis = seq_along(pos), eigenvalue = pos, pct_variance = 100 * pos / sum(pos))[1:10, ],
          file.path(OUT_TAB, "11_traitspace_axes.csv"))
site_scores <- bind_rows(lapply(c(Density = "den_mat", Biomass = "bio_mat"), function(nm) {
  W <- get(nm)[, tx_co$TAXANAME, drop = FALSE]; P <- W / rowSums(W)
  S <- P %*% as.matrix(tx_co[, c("A1", "A2")])
  tibble(SiteID = rownames(W), A1 = S[, 1], A2 = S[, 2])
}), .id = "Weighting")
write_csv(site_scores, file.path(OUT_TAB, "11_traitspace_site_scores.csv"))
say("  trait-space ordination written for step 06: %d taxa, axes 1 and 2 carry %.1f%% and %.1f%% of the variance",
    nrow(tx_co), 100 * pos[1] / sum(pos), 100 * pos[2] / sum(pos))

# ---- 4. Two trait groups that overlap: female dispersal and adult flying strength (24 Sep 2026) ---------------
# The chapter reports them as one result about dispersal capacity where both pass; the correlation of their codes
# across taxa is the record behind that sentence.
trf <- readRDS(IN$traits) %>% filter(TAXANAME %in% colnames(den_mat)) %>% as.data.frame()
pair <- function(a, b) { ok <- is.finite(trf[[a]]) & is.finite(trf[[b]]); tibble(state_a = a, state_b = b, n_taxa = sum(ok), pearson_r = cor(trf[[a]][ok], trf[[b]][ok]),
                                                                                     n_identical_codes = sum(abs(trf[[a]][ok] - trf[[b]][ok]) < 1e-12)) }
tgc <- bind_rows(pair("Disperse_High", "Flying_Strong"), pair("Disperse_Low", "Flying_Weak"))
write_csv(tgc, file.path(OUT_TAB, "11_trait_group_correlations.csv"))
say("  female dispersal against adult flying strength: high vs strong r = %.3f over %d taxa (%d identical codes)", tgc$pearson_r[1], tgc$n_taxa[1], tgc$n_identical_codes[1])

# ---- 5. Parent rows resolved before any metric (step 00c) ----------------------------------------------------
prov <- read_csv(file.path(RA, "tables", "H1", "H1_parent_rows_provenance.csv"), show_col_types = FALSE) %>% filter(mode == "resolution summary") %>% select(column, value = max_abs_diff)
above <- raw %>% group_by(TAXANAME) %>% summarise(sites = n_distinct(SiteID), .groups = "drop") %>%
  filter(grepl("idae$|inae$|ini$|oidea$|oidia$", TAXANAME) | TAXANAME %in% c("Acari","Hydracarina","Oligochaeta","Gastropoda","Turbellaria","Platyhelminthes","Collembola","Lepidoptera","Megaloptera","Thysanoptera","Chilopoda"))
prow <- bind_rows(prov, tibble(column = c("rows above genus rank still in the resolved table (parents with no child at those sites)", "sites at which such a row occurs"),
                               value = c(nrow(above), n_distinct(raw$SiteID[raw$TAXANAME %in% above$TAXANAME]))))   # 24 Sep night: distinct sites, not the sum of per-row site counts
write_csv(prow, file.path(OUT_TAB, "11_parent_rows_resolved.csv"))
say("  parent rows: %s parent-site cases resolved by DPAC at %s sites; %d rows above genus rank remain where no child was present",
    prov$value[prov$column == "parent-site cases resolved"], prov$value[prov$column == "sites with a resolution"], nrow(above))

# ---- 6. Presence-only trait space against richness, within region (24 Sep 2026) ---------------
# 28 Sep 2026: space-filling redundancy is S over FRic, so it is partly richness in another form; its correlation
# with richness, within region, is reported beside FRic's (supplement Table S4 part C; legacy Table S25).
fvr <- purrr::map_dfr(REGIONS, function(rg) { s <- dat %>% filter(Region == rg)
  s_fd <- metrics_raw$Richness_FD[match(s$SiteID, metrics_raw$SiteID)]   # the S in FR_SF (trait-set taxa present)
  tibble(Region = rg, n = nrow(s), rho_FRic_vs_richness = cor(s$FRic_den, s$Richness_alpha, method = "spearman"),
         rho_FRic_vs_trait_set_richness = cor(s$FRic_den, s_fd, method = "spearman"),
         rho_FRSF_vs_richness = cor(s$FR_SF_den_log10, s$Richness_alpha, method = "spearman"),
         rho_FRSF_vs_trait_set_richness = cor(s$FR_SF_den_log10, s_fd, method = "spearman"),
         rho_richness_vs_index = cor(s$BIBI, s$Richness_alpha, method = "spearman"), rho_FRic_vs_index = cor(s$BIBI, s$FRic_den, method = "spearman")) })
write_csv(fvr, file.path(OUT_TAB, "11_fric_vs_richness.csv"))
for (i in seq_len(nrow(fvr))) say("  %-15s FRic vs richness rho = %.3f (vs trait-set richness %.3f); FR_SF vs richness rho = %.3f (vs trait-set richness %.3f)",
                                  fvr$Region[i], fvr$rho_FRic_vs_richness[i], fvr$rho_FRic_vs_trait_set_richness[i],
                                  fvr$rho_FRSF_vs_richness[i], fvr$rho_FRSF_vs_trait_set_richness[i])

say("  trait-space fidelity (Maire et al. 2015): how well k axes reproduce Gower distance")
for (i in seq_len(nrow(dims))) {
  d <- dims[i, ]
  say("    %2d axes: %4.1f%% variance | mSD %.4f (Maire) | mAD %.4f (%.0f%% of mean distance) | r = %.3f",
      d$n_axes, d$cum_variance_pct, d$mSD_maire, d$mAD, d$mAD_pct_of_mean_distance, d$pearson_r)
}
r2 <- dims$pearson_r[dims$n_axes == 2]; r10 <- dims$pearson_r[dims$n_axes == 10]
say("  two axes reproduce trait dissimilarity at r = %.2f;", r2)
say("  ten axes at r = %.2f. The METRICS are unaffected -- dbFD computes", r10)
say("  them from the Gower distance or up to ten corrected axes, not from the two plotted axes,")
say("  so a two-axis ordination is a low-dimensional aid, not evidence.")
say("STEP 11 complete")
