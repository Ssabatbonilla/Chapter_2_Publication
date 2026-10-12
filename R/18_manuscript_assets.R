# =============================================================================
# 18_manuscript_assets.R
# =============================================================================
# WHAT THIS STEP DOES
#   Steps 02-17 write analysis output: raw CSVs at full numeric precision, named by the step that produced
#   them. This step assembles the manuscript set from them: formatted, numbered tables, the pipeline figures
#   copied under figure numbers, generated captions with every coefficient read from this run's tables, and a
#   manifest saying which claim each asset carries.
#
# INPUTS
#   output/tables/*.csv and output/figures/* from steps 01 to 20b; objects step 06 leaves in the same R session
#   (THRESH, HYP_ROWS, FIGS3_METRICS); the site table, land-use scores and H1_BIBI_thresholds_used.csv in the data
#   folder; Figure 1 (hand drawn) from data/hand_drawn_figures/ in the repository layout.
#
# NUMBERING
#   This step writes the LEGACY numbers (Tables 3a, 3b, 4, S1 to S29; Figures 1 to 6, S2, S3). Step 18b maps them
#   to the paper's numbers with R/supplement_asset_map.csv: Tables 3a and 3b become main-text Table 2, the kept
#   supplementary assets become Tables S1 to S8 and Figure S1 (legacy Figure S3), and Table 4, Tables S16, S17,
#   S19, S20, S28, Figure 6 and Figure S2 go to the companion trait-state manuscript. Figures 1 to 5 keep their
#   numbers.
#
# WHAT IT ENFORCES
#   - Every table rounded, with readable column names (fmt_table).
#   - Every figure drawn at journal final size (174 mm double column), 600 dpi,
#     Times New Roman where installed, 9.5 pt body, panel labels (a) (b) (c).
#   - A manifest mapping every asset to the claim it carries.
#
# OUTPUT
#   output/manuscript/figures/   Figure_<n>.png, .pdf and .svg (legacy numbers), as listed in ASSET_MANIFEST.csv,
#                                and Graphical_abstract.* when its source is found (see section 2b)
#   output/manuscript/tables/    Table_<n>.csv (legacy numbers), as listed in ASSET_MANIFEST.csv
#   output/manuscript/repository_tables/   verification and context tables, with INDEX.csv
#   output/manuscript/ASSET_MANIFEST.csv
#   output/manuscript/CAPTIONS.md
# =============================================================================

say("STEP 18 -- assemble manuscript-ready figures and tables")
suppressPackageStartupMessages(library(patchwork))

MS      <- file.path(OUT, "manuscript")
MS_FIG  <- file.path(MS, "figures"); MS_TAB <- file.path(MS, "tables")
for (d in c(MS, MS_FIG, MS_TAB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- 1. (removed) turnover versus nestedness ---------------------------------
# The beta-partitioning figure moved to its own paper on 8 September 2026. Step 16 still runs: the repository
# table regional_pools_shared_taxa.csv below reads 16_cross_region_pool_beta.csv.

# ---- 2. Table S1 (paper Table S1 part A) -- environmental and biological context by region ----
lu <- read_csv(IN$landuse, show_col_types = FALSE) %>%
  mutate(SiteID = std_site(SiteID)) %>% select(SiteID, ag_grad_lu)
cov <- read_csv(file.path(OUT_TAB, "10_site_trait_coverage.csv"), show_col_types = FALSE)
stopifnot(all(c("den_coverage", "bio_coverage") %in% names(cov)))

S1_CLASS <- { stopifnot(!anyNA(dat$BIBI_Class), all(dat$BIBI_Class %in% CLASS_LEVELS))
               tb <- table(dat$Region, factor(dat$BIBI_Class, levels = CLASS_LEVELS))
               z <- as.data.frame.matrix(tb); names(z) <- paste0("Streams rated ", tolower(CLASS_LEVELS))
               z$Region <- rownames(z); rownames(z) <- NULL; z[, c("Region", setdiff(names(z), "Region"))] }
T1 <- dat %>% left_join(lu, by = "SiteID") %>%
  left_join(cov %>% select(SiteID, den_coverage, bio_coverage), by = "SiteID") %>%
  group_by(Region) %>%
  summarise(`Sites (n)` = n(),
            `B-IBI median` = median(BIBI), `B-IBI min` = min(BIBI), `B-IBI max` = max(BIBI),
            `Selection gradient (land use), median` = median(ag_grad_lu, na.rm = TRUE),
            `Trait coverage, density (median %)` = 100 * median(den_coverage, na.rm = TRUE),
            `Trait coverage, biomass (median %)` = 100 * median(bio_coverage, na.rm = TRUE),
            `Total density median` = median(Total_Density),
            `Total biomass median` = median(Total_Biomass), .groups = "drop") %>%
  # 26 Sep 2026: the number of streams in each rating class, by region (regions are compared by condition class,
  # never by score), read from the class file step 01 checked against the release.
  left_join(S1_CLASS, by = "Region")
# 5 Oct 2026: the site totals are the median over the processed replicates of each
# replicate's total (the default), with the mean over the processed replicates beside them. When the site table carries
# the mean columns (R/00d_site_replicate_medians.R), Table S1 reports both bases, each as the regional median of the sites.
S1_SITE <- read_csv(IN$site_metrics, show_col_types = FALSE) %>% mutate(SiteID = std_site(SiteID))
if (all(c("Total_Density_mean", "Total_Biomass_mean") %in% names(S1_SITE))) {
  S1_MEAN <- dat %>% select(SiteID, Region) %>%
    left_join(S1_SITE %>% select(SiteID, Total_Density_mean, Total_Biomass_mean), by = "SiteID") %>%
    group_by(Region) %>%
    summarise(`Total density, mean of replicates (regional median, individuals m-2)` = median(Total_Density_mean),
              `Total biomass, mean of replicates (regional median, mg m-2)` = median(Total_Biomass_mean), .groups = "drop")
  T1 <- T1 %>%
    rename(`Total density, median of replicates (regional median, individuals m-2)` = `Total density median`,
           `Total biomass, median of replicates (regional median, mg m-2)` = `Total biomass median`) %>%
    left_join(S1_MEAN, by = "Region") %>%
    relocate(`Total density, mean of replicates (regional median, individuals m-2)`, .after = `Total density, median of replicates (regional median, individuals m-2)`) %>%
    relocate(`Total biomass, mean of replicates (regional median, mg m-2)`, .after = `Total biomass, median of replicates (regional median, mg m-2)`)
}
# 8 Oct 2026: the four rating-class cut-points of each region (every class is defined with its thresholds), read from
# the class file's threshold table (H1_BIBI_thresholds_used.csv), never typed. A region scored against two bioregions
# (Valley and Ridge: 29 streams SGV, one SRV) gets the bioregion that scores most of its streams in the four columns
# and the other in a text column, both from the same file.
S1_CUTS_RAW <- read_csv(file.path(H1, "H1_BIBI_thresholds_used.csv"), show_col_types = FALSE)
stopifnot(all(c("Typology", "Label", "Score", "Bioregion", "n_sites") %in% names(S1_CUTS_RAW)),
          setequal(unique(S1_CUTS_RAW$Typology), REGIONS),
          all(table(S1_CUTS_RAW$Typology, S1_CUTS_RAW$Bioregion)[table(S1_CUTS_RAW$Typology, S1_CUTS_RAW$Bioregion) > 0] == 4))
S1_CUT_LABELS <- c("VeryPoor–Poor" = "Cut-point, very poor | poor", "Poor–Fair" = "Cut-point, poor | fair",
                   "Fair–Good" = "Cut-point, fair | good", "Good–Excellent" = "Cut-point, good | excellent")
stopifnot(setequal(unique(S1_CUTS_RAW$Label), names(S1_CUT_LABELS)))
S1_BR <- S1_CUTS_RAW %>% distinct(Typology, Bioregion, n_sites) %>% group_by(Typology) %>%
  arrange(desc(n_sites), .by_group = TRUE) %>% mutate(main = row_number() == 1) %>% ungroup()
stopifnot(all((S1_BR %>% group_by(Typology) %>% summarise(n = sum(n_sites)))$n == 30))
S1_MAIN <- S1_CUTS_RAW %>% semi_join(S1_BR %>% filter(main), by = c("Typology", "Bioregion")) %>%
  mutate(Label = unname(S1_CUT_LABELS[Label])) %>%
  tidyr::pivot_wider(id_cols = c(Typology, Bioregion, n_sites), names_from = Label, values_from = Score) %>%
  transmute(Region = Typology, `Bioregion of the cut-points (streams scored)` = sprintf("%s (%d)", toupper(Bioregion), n_sites),
            !!!rlang::syms(unname(S1_CUT_LABELS)))
S1_OTHER <- S1_CUTS_RAW %>% semi_join(S1_BR %>% filter(!main), by = c("Typology", "Bioregion")) %>%
  mutate(Label = factor(Label, levels = names(S1_CUT_LABELS))) %>% arrange(Typology, Bioregion, Label) %>%
  group_by(Region = Typology, Bioregion, n_sites) %>%
  summarise(txt = sprintf("%s (%d stream%s): %s", toupper(first(Bioregion)), first(n_sites), ifelse(first(n_sites) == 1, "", "s"),
                          paste(sprintf("%.2f", Score), collapse = ", ")), .groups = "drop") %>%
  group_by(Region) %>% summarise(`Other bioregion and its cut-points` = paste(txt, collapse = "; "), .groups = "drop")
T1 <- T1 %>% left_join(S1_MAIN, by = "Region") %>% left_join(S1_OTHER, by = "Region") %>%
  mutate(`Other bioregion and its cut-points` = coalesce(`Other bioregion and its cut-points`, "none"))
stopifnot(!anyNA(T1[["Cut-point, poor | fair"]]))
say("  Table S1 cut-points: %s", paste(sprintf("%s %s", T1$Region, T1[["Bioregion of the cut-points (streams scored)"]]), collapse = "; "))
T1 <- T1 %>% fmt_table(2)

# ---- 3. Assemble every table, formatted -------------------------------------
rd <- function(f) read_csv(file.path(OUT_TAB, f), show_col_types = FALSE)

# The metric dictionary, emitted so the reader sees the same definitions the
# figures are labeled from: name, meaning, equation where one is needed, and the source.
# "Computed on" says, per index, which distance or space the code used (24 Sep 2026; checked against
# R/13_fd_provenance.R, R/00b_redundancy_geometry.R and the FD::dbFD source). No single corrected matrix serves
# every metric: RaoQ, FDis and overlap redundancy use the uncorrected Gower distance.
GOWER <- "group-weighted Gower distance over the 41 states (each trait group weighted by 1/its number of states; a missing group dropped for the pair; FD::gowdis)"
COMPUTED_ON <- c(
  FRic_den = paste0("Convex hull volume in the first 10 axes of a Cailliez-corrected principal coordinate analysis of the ", GOWER,
                    sprintf("; standardized by the volume of all %d taxa in the trait set (FD::dbFD)", ncol(den_mat))),
  FEve_den = "Minimum spanning tree on all axes of the Cailliez-corrected principal coordinate analysis (FD::dbFD)",
  FEve_bio = "Minimum spanning tree on all axes of the Cailliez-corrected principal coordinate analysis (FD::dbFD)",
  FDiv_den = "First 10 axes of the Cailliez-corrected principal coordinate analysis (FD::dbFD)",
  FDiv_bio = "First 10 axes of the Cailliez-corrected principal coordinate analysis (FD::dbFD)",
  FDis_den = paste("The uncorrected", GOWER, "(FD::dbFD)"),
  FDis_bio = paste("The uncorrected", GOWER, "(FD::dbFD)"),
  RaoQ_den = paste0("The uncorrected ", GOWER, "; half the weighted mean squared distance (FD::dbFD)"),
  RaoQ_bio = paste0("The uncorrected ", GOWER, "; half the weighted mean squared distance (FD::dbFD)"),
  FR_RD_den = "1 - Q/D, with Q the weighted mean unsquared Gower distance (not the RaoQ above) and D Simpson diversity (R/00b)",
  FR_RD_bio = "1 - Q/D, with Q the weighted mean unsquared Gower distance (not the RaoQ above) and D Simpson diversity (R/00b)",
  FR_SF_den_log10 = "log10(S / FRic + 1), S the number of trait-set taxa present (R/00b)",
  FR_DB_RD = "Overlap redundancy weighted by biomass minus overlap redundancy weighted by density (R/00b)")
T_DICT <- METRIC_DICT %>%
  mutate(computed_on = unname(COMPUTED_ON[key]), computed_on = ifelse(is.na(computed_on), "", computed_on)) %>%
  transmute(Property = unname(c(structure = "Structure and tolerance", totals = "Totals", `functional diversity` = "Functional diversity",
                                redundancy = "Redundancy", inequality = "Inequality")[property]),   # one display name per family
            Metric = name, Symbol = symbol, Weighting = currency,
            `What it measures` = meaning, Equation = equation, `Computed on` = computed_on, Source = citation) %>%
  distinct(Property, Metric, Symbol, Weighting, .keep_all = TRUE)
stopifnot(all(names(COMPUTED_ON) %in% METRIC_DICT$key))

# Legacy asset numbers; step 18b maps them to the paper's (R/supplement_asset_map.csv).
# Authored tables with no pipeline source are not written here: legacy Table S2 (sampling design, paper Table S1
# part B) and Table S3 (trait groups, paper Table S2 part A). Beta partitioning is not a Chapter 2 asset.
# Tables 3a, 3b and 4 are derived here rather than read whole: 3a and 3b split the responsive-metric screen by
# description (step 18b merges them into the paper's Table 2), and Table 4 is the per-block fourth-corner record
# (companion trait-state manuscript).
T2_ALL     <- rd("07_T2_responsive_metrics.csv")
# Standing stock belongs to Table 3b, not Table 3a: Table 3a holds the structural metrics only and the two totals
# sit in Table 3b (9 September 2026). The filter reads the family display name step 07 writes (its `property` map;
# "Structure and tolerance" since 25 Sep 2026), and the step stops if Table 3a is empty or the label is absent.
STRUCTURAL <- c("Structure and tolerance")
T3A <- dplyr::filter(T2_ALL,  Property %in% STRUCTURAL)
T3B <- dplyr::filter(T2_ALL, !Property %in% STRUCTURAL)
stopifnot("Table 3a has no rows: the structural family label has changed" = nrow(T3A) > 0,
          all(STRUCTURAL %in% unique(T2_ALL$Property)))

# The per-block fourth-corner record (legacy Table 4; companion trait-state manuscript), carried under both
# weightings, so Currency is a column and not a file.
T4 <- dplyr::bind_rows(
        dplyr::mutate(rd("15_fourthcorner_per_block.csv"),         Currency = "Density"),
        dplyr::mutate(rd("15_fourthcorner_per_block_biomass.csv"), Currency = "Biomass")) %>%
      # A state passes only under every one of the five seeds (22 Sep 2026); a state that
      # passes under some seeds is listed as borderline with its FDR range, never as passing.
      dplyr::filter(n_seeds_pass > 0) %>%
      dplyr::mutate(`FDR p, range over seeds` = sprintf("%.3f to %.3f", p_fdr_min, p_fdr_max),
                    `Seeds passing` = sprintf("%d of %d", n_seeds_pass, N_SEEDS_FC)) %>%
      dplyr::arrange(dplyr::desc(n_seeds_pass == N_SEEDS_FC), Currency, Region, block) %>%
      dplyr::select(Region, Weighting = Currency, `Trait group` = block, `Trait state` = term,
                    Statistic = stat, `raw p` = p_raw, `FDR p` = p_fdr, `FDR p, range over seeds`,
                    `Seeds passing`, `Carrier taxa` = n_carriers, Verdict = verdict)

TABLES <- list(
  list("Table_3a", "Structural metrics meeting the screening criterion, by region", T3A, "H1"),
  list("Table_3b", "Trait, total, inequality and weighting metrics meeting the screening criterion, by region",
       T3B, "H2, H3"),
  # The record behind legacy Figure 6 (fig_traitstate_tiers); statistics in legacy Table S16. All three go to the
  # trait-state manuscript (step 18b).
  list("Table_4",  "Trait states associated with condition, by trait group, region and weighting, corrected within each region-by-trait-group block (each permutation model adjusted within the block and the larger adjusted p reported, ade4 model 6; passing under all five seeds; borderline states marked; carrier taxa per state); the record behind Figure 4",
       T4, "H3, mechanism"),
  # Legacy Table_S2 (sampling design; paper Table S1 part B) and Table_S3 (trait groups and states; paper Table S2
  # part A) are authored in the supplement and have no pipeline source.
  list("Table_S1", "Environmental and biological context by region, with the index's rating-class cut-points. Cut-points are each bioregion's 50th, 25th and 10th percentiles of reference-site scores and half the 10th (Smith et al., 2017, Table 5); a score equal to a cut-point takes the better class. SGV is the Southern Great Valley metric set and SRV the Southern Ridge and Valley set; UNP and MAC are the bioregions that score the Piedmont and Coastal Plain streams", T1, "Study system; Methods 2.2 (rating classes)"),
  # 28 Sep 2026: one coverage table (paper Table S3 part A). The 80% subset (step 10) and the partial correlation
  # holding coverage constant (step 12) joined by metric and region; legacy S5 keeps its source column names. Legacy S4,
  # S6, S8, S12, S13 and S22 are written to repository_tables/ below, under descriptive names.
  list("Table_S5", "Trait coverage and the responsive associations: each association in the streams meeting the 80% trait-coverage criterion, and partialled on the share of density or biomass in taxa with a trait record",
       dplyr::full_join(rd("10_coverage80_sensitivity.csv"),
                        rd("12_coverage_confound_test.csv") %>% dplyr::rename(Metric_family = Family),
                        by = c("Metric", "Region")), "limitation (coverage; T8)"),
  list("Table_S7", "Metric definitions, weighting and sources", T_DICT, "Methods 2.4"),
  # Class-based within-region tests (PERMANOVA, ANOSIM and dispersion on condition class) are dropped (24 Sep
  # 2026): condition enters as the continuous index and classes label figures only, and two regions have classes of
  # one site. Step 04 still writes them to its own output.
  list("Table_S9", "Multivariate tests, cross-region and within-region",
       rd("07_S4_multivariate_tests.csv") %>% dplyr::filter(!(Scope == "within-region" & Term == "B-IBI class")), "C6"),
  list("Table_S10","Full metric screen, all metrics and regions", rd("07_S1_full_screen.csv"), "C1"),
  # Legacy Table_S11 (land-use covariate models) was removed 24 Sep 2026: no Chapter 2 analysis uses land use.
  # The legacy numbers are kept here; step 18b assigns the paper's numbers.
  # The Sorensen partition columns are dropped from the regional-pool table (24 Sep 2026): beta diversity moved to its
  # own paper. Step 16 still writes them to its own output.
  list("Table_S14","Minimum detectable effect and achieved power",
       rd("14_achieved_power_by_metric.csv"), "C2"),
  list("Table_S15", sprintf("Null-model standardized effect sizes (trait labels shuffled, %d permutations)", N_PERM_SES), rd("17_ses_vs_raw.csv"), "C3"),
  list("Table_S16","Fourth-corner results by trait group, both weightings, corrected within each region-by-trait-group block (each permutation model adjusted within the block and the larger adjusted p reported, ade4 model 6), with the correction across all states of a region (Benjamini-Hochberg on the model-6 p) beside each state",
       dplyr::bind_rows(dplyr::mutate(rd("15_fourthcorner_per_block.csv"), Weighting = "Density"),
                        dplyr::mutate(rd("15_fourthcorner_per_block_biomass.csv"), Weighting = "Biomass")) %>%
         dplyr::left_join(dplyr::bind_rows(dplyr::mutate(rd("15_fourthcorner_regionwide.csv"), Weighting = "Density"),
                                           dplyr::mutate(rd("15_fourthcorner_regionwide_biomass.csv"), Weighting = "Biomass")) %>%
                            dplyr::select(Weighting, Region, term, n_tests_regionwide, n_seeds_pass_regionwide, p_fdr_regionwide_min, p_fdr_regionwide_max, verdict_regionwide),
                          by = c("Weighting", "Region", "term")) %>%
         dplyr::relocate(Weighting), "C8"),
  # Added 24 September 2026: the multivariate per-group fourth-corner test, the leverage audit, leave-one-taxon-out
  # and the density-against-biomass difference test. Legacy numbers; step 18b assigns the paper's.
  list("Table_S17","One multivariate fourth-corner test per trait group (RLQ inertia under permutation model 6), corrected across the trait groups of each region, with the group's states in a fuzzy correspondence analysis and with the states standardized",
       dplyr::bind_rows(dplyr::mutate(rd("15_fourthcorner_pergroup.csv"), Weighting = "Density"),
                        dplyr::mutate(rd("15_fourthcorner_pergroup_biomass.csv"), Weighting = "Biomass")) %>% dplyr::relocate(Weighting), "L37 sensitivity"),
  list("Table_S18","Leave-few-out leverage audit of every responsive screen association", rd("03b_leverage_audit.csv"), "Methods 2.5 (R8)"),
  list("Table_S19","Leave-one-taxon-out for every trait state passing under all five seeds",
       dplyr::bind_rows(dplyr::mutate(rd("15_leave_one_taxon_out.csv"), Weighting = "Density"),
                        dplyr::mutate(rd("15_leave_one_taxon_out_biomass.csv"), Weighting = "Biomass")) %>% dplyr::relocate(Weighting), "FL1"),
  list("Table_S20","Density-against-biomass difference of each trait state's fourth-corner statistic, streams and taxa resampled independently with replacement on each draw",
       rd("15b_weighting_difference.csv"), "H3"),
  # 26 Sep 2026: the taxon-level screen and the tolerance of the responding taxa (paper Table S8 parts A and B).
  list("Table_S23","Taxon-level screen: every taxon occurring at eight or more sites in a region against the within-region B-IBI gradient, weighted by density and by biomass",
       rd("09_taxon_response_screen.csv") %>% dplyr::rename(Weighting = Currency, Taxon = TAXANAME, `Sites (n)` = n_occ,
                                                           `Mean share` = mean_share, `raw p` = p_raw, `FDR p` = p_FDR), "Results 3.2 (H1)"),
  list("Table_S24","Tolerance values of the taxa whose density or biomass was associated with condition",
       rd("09b_taxon_tolerance_screen.csv") %>% dplyr::filter(Response %in% c("increaser", "decreaser")) %>%
         dplyr::select(Region, Weighting = Currency, Taxon = TAXANAME, Family, Response, rho, `FDR p` = p_FDR,
                       `Tolerance value` = TV, `Tolerance class` = TV_class), "Results 3.2 (H1)"),
  # 28 Sep 2026: the presence-only controls against richness, within region, so the reader can see how much of
  # functional richness and of space-filling redundancy (S / FRic) is richness itself (paper Table S4 part C).
  list("Table_S25","Functional richness and space-filling redundancy against taxonomic richness and trait-set richness within each region (Spearman correlation across streams)",
       rd("11_fric_vs_richness.csv") %>%
         dplyr::transmute(Region, `Streams (n)` = n,
                          `FRic vs taxonomic richness` = rho_FRic_vs_richness, `FRic vs trait-set richness` = rho_FRic_vs_trait_set_richness,
                          `FR_SF vs taxonomic richness` = rho_FRSF_vs_richness, `FR_SF vs trait-set richness` = rho_FRSF_vs_trait_set_richness),
       "Methods 2.4; Appendix S1 (the weighting-invariant controls)"),
  # 8 Oct 2026 (step 20): the sensitivity of the density-weighted trait results to the four chironomid groups
  # (paper Table S7 parts A to C; legacy Table S28 goes to the trait-state manuscript).
  list("Table_S26","Share of the trait-set total made up by the four chironomid taxa (the tribes Chironomini and Tanytarsini and the subfamilies Orthocladiinae and Tanypodinae) and its association with condition, by region and weighting",
       rd("20_chironomid_share.csv"), "Results 3.3 (chironomid sensitivity); Methods 2.5"),
  list("Table_S27","Weighted functional metrics and the redundancy divergence with every taxon and without the four chironomid taxa, against removing the four most abundant other taxa and 1,000 random sets of four other taxa (Spearman correlation with B-IBI, within region; false discovery rate within metric family and region as in Table S10)",
       rd("20_chironomid_functional_metrics.csv") %>%
         dplyr::transmute(Weighting, Metric, Region, `rho, every taxon` = rho_all, `FDR p, every taxon` = fdr_all, `Meets the criterion, every taxon` = met_all,
                          `rho, without the four` = rho_removed, `FDR p, without the four` = fdr_removed, `Meets the criterion, without the four` = met_removed,
                          `rho, four most abundant other taxa removed` = rho_top4_other, `Those taxa` = top4_other,
                          `rho, random sets of four, 5th percentile` = random4_p05, `rho, random sets of four, 95th percentile` = random4_p95,
                          `Random sets moving rho as far as removing the four` = random4_as_extreme),
       "Results 3.3 (chironomid sensitivity); Methods 2.5"),
  list("Table_S28","Trait states passing under all five seeds with every taxon or without the four chironomid taxa (fourth-corner, permutation model 6, 9,999 permutations, false discovery rate within region and trait group)",
       rd("20_chironomid_fourthcorner_states.csv") %>%
         dplyr::transmute(Weighting, Region, `Trait group` = unname(TRAIT_GROUP_SHORT[block]), `Trait state` = unname(TRAIT_STATE_LAB[state]),
                          `Statistic, every taxon` = stat_kept, `FDR p (largest of five seeds), every taxon` = p_fdr_max_kept, `Verdict, every taxon` = verdict_all,
                          `Statistic, without the four` = stat_removed, `FDR p (largest of five seeds), without the four` = p_fdr_max_removed, `Verdict, without the four` = verdict_removed,
                          `Share of the state's weighted mass made up by the four (%)` = 100 * chiro_part_of_mass),
       "Results 3.3 (chironomid sensitivity); Methods 2.5"),
  list("Table_S29","Genus-level split of the four chironomid taxa: Rao's quadratic entropy and overlap redundancy weighted by density against condition under each scenario (1,000 draws per scenario and region)",
       rd("20_chironomid_genus_mix_simulation.csv") %>%
         dplyr::transmute(Region, Metric, Scenario, `rho, median` = rho_median, `rho, 5th percentile` = rho_p05, `rho, 95th percentile` = rho_p95,
                          `Draws with rho of at least 0.30 in the direction of the group-profile analysis` = share_in_baseline_direction_abs_ge_030),
       "Methods 2.5; Discussion 4.5 and 4.6")
  # The taxon-level response screen is legacy Table_S23 above (written by step 09).
)

# 28 Sep 2026: verification, provenance and context tables leave the supplement for the repository. They are
# still written, under descriptive names, so the text's numbers keep a source; INDEX.csv says what each was.
REPO_TAB <- file.path(MS, "repository_tables"); dir.create(REPO_TAB, showWarnings = FALSE)
REPO_TABLES <- list(
  list("trait_coverage_by_site", "Trait coverage by site and region", rd("10_site_trait_coverage.csv"), "Table S4"),
  list("traitspace_dimensionality", "Trait-space dimensionality (Maire et al. 2015)", rd("11_traitspace_dimensionality.csv"), "Table S6"),
  list("fd_metric_provenance", "Functional-diversity metric provenance verification", rd("13_fd_provenance_verification.csv"), "Table S8"),
  list("regional_pools_shared_taxa", "Cross-region pool comparison: taxa per region and shared taxa",
       rd("16_cross_region_pool_beta.csv") %>% dplyr::select(-dplyr::any_of(c("turnover", "nestedness", "total", "pct_nestedness", "turnover_exceeds_nestedness"))), "Table S12"),
  list("regional_metric_distributions", "Regional metric distributions", rd("07_S2_regional_distributions.csv"), "Table S13"),
  list("parent_rows_resolved", "Parent rows resolved before any metric, distributed among child taxa in proportion to their density (DPAC)",
       read_csv(file.path(RA, "tables", "H1", "H1_parent_rows_resolution.csv"), show_col_types = FALSE), "Table S22"),
  # 8 Oct 2026: the record behind Discussion 4.5's "Chironomidae made up the largest share of density in every region".
  list("family_density_biomass_shares", "Each family's share of regional density and biomass, every taxon, with the number of streams in which it has the largest share of density",
       rd("10_family_density_share.csv"), "new, 8 Oct 2026 (Discussion 4.5)")
)
repo_idx <- purrr::map_dfr(REPO_TABLES, function(z) {
  df <- fmt_table(z[[3]]); write_csv(df, file.path(REPO_TAB, paste0(z[[1]], ".csv")))
  tibble(File = paste0(z[[1]], ".csv"), Title = z[[2]], Formerly = z[[4]], Rows = nrow(df), Cols = ncol(df)) })
write_csv(repo_idx, file.path(REPO_TAB, "INDEX.csv"))
say("  repository tables: %d written to %s (formerly %s)", nrow(repo_idx), if (REPO_LAYOUT) "output/manuscript/repository_tables" else REPO_TAB, paste(repo_idx$Formerly, collapse = ", "))

tab_rows <- purrr::map_dfr(TABLES, function(z) {
  df <- fmt_table(z[[3]])
  write_csv(df, file.path(MS_TAB, paste0(z[[1]], ".csv")))
  tibble(Asset = z[[1]], Kind = "table",
         Placement = ifelse(grepl("_S", z[[1]]), "supplementary", "main"),
         Title = z[[2]], Carries = z[[4]], Rows = nrow(df), Cols = ncol(df))
})

# ---- 4. Assemble figures with final numbering -------------------------------
# Sources are the pipeline figure stems; the Figure_N names are the legacy manuscript numbering (step 18b maps the
# supplementary ones). Figure 2 is the study area map from step 19.
#
# Figure 1 is drawn by hand and has no pipeline source, so it is copied in rather than built: from
#   data/hand_drawn_figures in the public repository, otherwise from the author's project tree. It is exported at
#   the specification the pipeline enforces for everything else (174 mm, 600 dpi, Times New Roman).
#
# Not carried: fig4_traitspace_binned, figS3_structure_scatters, figS4_function_scatters and the working figures
# of step 06.
FIG1_SRC <- if (REPO_LAYOUT) file.path(DATA, "hand_drawn_figures", "Figure_1") else
            file.path(CH2_ROOT, "Manuscripts", "Figure1_Conceptual", "Figure_1")

FIGS <- list(
  list("Figure_1", NA_character_, "main",
       "Three ways a stream's macroinvertebrates respond along the biological condition gradient, drawn as hypotheses 1 to 3",
       "H1, H2, H3"),
  list("Figure_2", "fig1_study_area", "main",
       "Study area and the condition gradient each region spans", "design"),
  # Figures 3 (hypotheses 1 and 2) and 4 (hypothesis 3) are the responsive metrics by hypothesis; Figure 5 is the
  # effect-size plot (a, b).
  list("Figure_3", "fig3_hyp12_responsive", "main",
       "Every metric of hypotheses 1 and 2 that met the criterion in at least one region, against the condition gradient within region, both weightings overlaid, with each region's Poor | Fair cut", "H1, H2"),
  list("Figure_4", "fig4_hyp3_responsive", "main",
       "Every metric of hypothesis 3 that met the criterion in at least one region: the two weightings, the totals, inequality and the two divergences, within region", "H3"),
  list("Figure_5", "fig5_effect_sizes_by_region", "main",
       "Every metric-by-weighting test against the condition gradient, within region", "H1, H2, H3"),
  # Legacy Figure 6 and Figure S2: the trait-state figures, moved by step 18b to the companion trait-state manuscript.
  list("Figure_6", "fig_traitstate_tiers", "main",
       "Trait states associated with condition within their trait groups, and the stricter corrections each also met, by region and weighting", "H3, mechanism"),
  # The binned trait-space figure (fig4_traitspace_binned) is not carried: it carries no claim.
  list("Figure_S2", "figS4_traitstate_difference", "supplementary",
       "The density-against-biomass difference of each trait state's fourth-corner statistic, by region", "H3, record"),
  # Steps 20 and 20b; paper Figure S1 after step 18b.
  list("Figure_S3", "figS5_chironomid_removal", "supplementary",
       "Density-weighted functional metrics and trait states with every taxon and with the four chironomid groups removed", "Results 3.3 (chironomid sensitivity)")
)

fig_rows <- purrr::map_dfr(FIGS, function(z) {
  ok <- FALSE
  # Figure 1 is hand-drawn: its source is the exported conceptual figure, not a
  # pipeline stem.
  stem <- if (is.na(z[[2]])) FIG1_SRC else file.path(OUT_FIG, z[[2]])
  unlink(file.path(MS_FIG, paste0(z[[1]], c(".png", ".pdf", ".svg"))))     # no stale copy survives a source that lost a format
  for (ext in c(".png", ".pdf", ".svg")) {
    src <- paste0(stem, ext)
    if (file.exists(src)) { file.copy(src, file.path(MS_FIG, paste0(z[[1]], ext)), overwrite = TRUE); ok <- TRUE }
  }
  tibble(Asset = z[[1]], Kind = "figure", Placement = z[[3]],
         Title = z[[4]], Carries = z[[5]], Rows = NA_integer_, Cols = NA_integer_,
         built = ok)
})

# ---- 2b. Graphical abstract (6 Oct 2026) ---------------------------------------------------------------------------
# The journal graphical abstract (slope pairs by region, 80 x 80 mm, 600 dpi) is drawn by a separate script from this
# pipeline's tables, outside this repository. It is not a numbered figure. When its files are found at GA_SRC they are
# copied here as Graphical_abstract.{svg,pdf,png} so the submission package holds them, and its caption is generated
# below with every coefficient read from the screen table. GA_SRC is a path in the author's project tree; when it is
# absent the asset is skipped with a notice.
GA_SRC <- if (REPO_LAYOUT) NA_character_ else   # the repository never reads outside itself; its copy of the asset ships in output/
  file.path(CH2_ROOT, "Manuscripts", "Graphical_abstract_2026-10-06", "stage5", "candidates", "J_B_slope_pairs")
if (!REPO_LAYOUT) unlink(file.path(MS_FIG, paste0("Graphical_abstract", c(".png", ".pdf", ".svg"))))   # the repository cannot redraw it, so it keeps the shipped copy
ga_ok <- FALSE
for (ext in c(".png", ".pdf", ".svg")) {
  src <- paste0(GA_SRC, ext)
  if (file.exists(src)) { file.copy(src, file.path(MS_FIG, paste0("Graphical_abstract", ext)), overwrite = TRUE); ga_ok <- TRUE }
}
if (!ga_ok) say("  NOTICE: graphical abstract %s; the asset is not written", if (REPO_LAYOUT) "is drawn outside the repository (its copy ships in output/manuscript/figures/)" else paste("source not found:", GA_SRC))
fig_rows <- bind_rows(fig_rows, tibble(
  Asset = "Graphical_abstract", Kind = "figure", Placement = "graphical abstract",
  Title = "Trait sharing weighted by density and by biomass in the poorer and the better third of each region's condition gradient (J-B, slope pairs by region)",
  Carries = "H3", Rows = NA_integer_, Cols = NA_integer_, built = ga_ok))

# An asset the manuscript numbers but the asset set does not hold is a hole in
# the submission package, and an asset the set holds but the manuscript does not
# number is an orphan. Report both; neither is allowed to pass silently.
#
# TABLES ARE CHECKED TOO, and they are the more dangerous half. A stale figure
# file sits under a number the manuscript does not use, so it is merely extra.
# A stale TABLE file can sit under a number the manuscript DOES use, for
# different content, after a renumbering.
check_orphans <- function(dir, expected, what) {
  found <- unique(sub("[.][^.]+$", "", list.files(dir)))
  orphans <- setdiff(found, expected)
  if (length(orphans))
    say("  WARNING: %d orphaned %s file(s) in %s: %s", length(orphans), what,
        basename(dir), paste(sort(orphans), collapse = ", "))
  holes <- setdiff(expected, found)
  if (length(holes))
    say("  WARNING: %d %s(s) the manuscript numbers but this set does not hold: %s",
        length(holes), what, paste(sort(holes), collapse = ", "))
  if (!length(orphans) && !length(holes))
    say("  %s set complete: %d assets, no orphans, no holes", what, length(expected))
}
check_orphans(MS_FIG, c(vapply(FIGS,   `[[`, character(1), 1), "Graphical_abstract"), "figure")
check_orphans(MS_TAB, vapply(TABLES, `[[`, character(1), 1), "table")

missing <- fig_rows %>% filter(!built)
if (nrow(missing))
  say("  WARNING: %d figure(s) not found in output/figures: %s",
      nrow(missing), paste(missing$Asset, collapse = ", "))

manifest <- bind_rows(tab_rows, fig_rows %>% select(-built)) %>%
  arrange(Placement, Kind, as.integer(gsub("\\D", "", Asset)), Asset)     # numeric order: S4 before S10
# Table 4 is the record behind legacy Figure 6, not a main table.
manifest$Placement[manifest$Asset == "Table_4"] <- "record"
# The Caption column is filled in and the manifest written after section 5, once
# it is known which captions were actually found in the manuscript.

# ---- 5. Captions, generated from the same tables the figures are drawn from --
#
# WHY THESE ARE GENERATED, NOT WRITTEN BY HAND
#   Statistics were removed from the figure panels, where they overlapped the
#   data and could drift out of step with the screen. They have to reappear
#   somewhere, and the caption is the right place: it is where a reader looks
#   for them, it costs no space in the graphic, and it can be edited without
#   redrawing anything. Generating it from `screen` means a caption can never
#   disagree with the figure it describes -- the same guarantee step 06 gives
#   for the figures themselves.
#
# THE CAPTION PATTERN
#   1. What the figure shows, and the screening criterion.
#   2. Per panel: label, metric abbreviation, every responsive coefficient, and
#      ONE interpretation clause saying what the direction means.
#   3. What the symbols, colors and lines encode.
#   4. Where the rest of the evidence is.
#   Interpretation clauses avoid spatial-motion constructions, and the study
#   organisms are macroinvertebrates.

stat_str <- function(metric) {
  r <- screen %>% filter(Metric == metric, Tier == "strong") %>% arrange(Region)
  if (!nrow(r)) return("no region met the screening criterion")
  paste(sprintf("%s rho = %+.2f, FDR %s", r$Region, r$rho,
                ifelse(r$p_FDR < 0.001, "< 0.001", sprintf("= %.3f", r$p_FDR))),
        collapse = "; ")
}


panel_caption <- function(letter, label, metric, meaning) {
  # Label and abbreviation both come from METRIC_DICT so a caption cannot use a
  # different name from the axis it describes.
  label <- mdict(metric, "name")
  ab    <- mdict(metric, "symbol")
  if (identical(ab, "-")) ab <- label
  # Do not print "total density (total density; ...)" -- where the abbreviation
  # IS the label, give the statistics alone.
  head <- if (identical(tolower(ab), tolower(label)))
            sprintf("(%s) %s (%s)", letter, label, stat_str(metric))
          else sprintf("(%s) %s (%s; %s)", letter, label, ab, stat_str(metric))
  # Guarantee terminal punctuation so the joined list never runs into the
  # sentence that follows it.
  if (!grepl("[.]$", meaning)) meaning <- paste0(meaning, ".")
  paste(head, "--", meaning)
}

CRITERION <- sprintf(
  "Metrics were screened within region against the B-IBI gradient using Spearman rank correlation, with a metric counted as responsive at |rho| >= %.2f and Benjamini-Hochberg adjusted p <= %.2f.",
  RHO_MIN, FDR_MAX)
# Axis context, generated from Y_DOMAIN so the caption cannot describe an axis
# the figure does not have. This is the answer to the two questions a reader
# asks of every one of these panels: what range is possible, and why does the
# axis not start at zero.
axis_note <- function(keys) {
  syms <- unique(vapply(keys, function(k) mdict(k, "symbol"), character(1)))
  syms <- syms[syms %in% names(Y_DOMAIN)]
  if (!length(syms)) return("")
  parts <- vapply(syms, function(sy) {
    d <- Y_DOMAIN[[sy]]
    sprintf("%s %s to %s", sy, format(d[1], trim = TRUE), format(d[2], trim = TRUE))
  }, character(1))
  paste0("Each metric is drawn on a fixed axis, identical across regions and ",
         "across density and biomass weighting, so panels of the same metric ",
         "are directly comparable: ", paste(parts, collapse = "; "), ". ",
         "Axes are set to round bounds that bracket the full observed range ",
         "rather than to each panel's extremes, and no observation falls ",
         "outside them.")
}

BOUNDED_NOTE <- paste(
  "Functional divergence, overlap redundancy and the Gini coefficient are",
  "bounded on [0, 1] by construction, and the redundancy and concentration",
  "divergences on [-1, 1]. Their axes start above the theoretical minimum",
  "because the observed values occupy a narrow part of that range and a full",
  "0-1 axis would compress the within-region gradient to invisibility;",
  "Rao's quadratic entropy, functional dispersion, richness and the",
  "Hilsenhoff Biotic Index, which have a meaningful zero and data reaching",
  "toward it, are drawn from zero.")

SYMBOLS <- "Point color gives physiographic region and point shape gives B-IBI condition class. Fitted lines are ordinary least-squares fits with 95% confidence bands, drawn only for regions in which the metric met the screening criterion, so a line is never shown where the screen did not support an association."

# cap2 (the caption of the structure scatter grid, figS3_structure_scatters) is built from FIGS3_METRICS (step 06)
# but is not written: that figure is not carried.
stopifnot(exists("FIGS3_METRICS"))
# 3 Oct 2026. The threshold clause reads the cut from the same table step 06 draws it from.
stopifnot(exists("THRESH"))
THRESH_NOTE <- sprintf(
  "The vertical line in each panel is the region's Poor | Fair rating threshold at the family level (Smith et al., 2017; %s), drawn as a label of the index's own classes and not as a test; streams to its left are rated very poor or poor, to its right fair or better.",
  paste(sprintf("%s %.2f", THRESH$Region, THRESH$cut), collapse = ", "))
VR_TWO_BIOREGIONS <- "Valley and Ridge is scored against two bioregions, and the line is the one that scores 29 of its 30 streams (Data/Data/build_bibi_site_classes.R)."
cap2 <- paste(
  "**Figure S3.** Structure and totals along the within-region B-IBI gradient.",
  CRITERION,
  sprintf("Rows are %s; columns are regions; each point is one stream, and the y axis of each row is its own.",
          paste(vapply(FIGS3_METRICS, function(m) mdict(m, "name"), character(1)), collapse = ", ")),   # names as the dictionary spells them (Shannon, Hilsenhoff are proper nouns; 3 Oct 2026)
  SYMBOLS,
  THRESH_NOTE,
  "Full screening results for all metrics and regions are in Table S10.")

# values the generated captions of legacy Figures 6 and S2 quote
WD_PER_REGION <- { z <- rd("15b_weighting_difference.csv") %>% dplyr::filter(keep) %>% dplyr::count(Region); stopifnot(length(unique(z$n)) == 1); unique(z$n) }
WD_RESAMPLES  <- as.integer(rd("15b_weighting_difference_settings.csv") %>% dplyr::filter(setting == "resamples per seed") %>% dplyr::pull(value))
FC_FDR        <- format(FDR_MAX, nsmall = 2)
N_RW_STATES   <- { z <- unique(rd("15_fourthcorner_regionwide.csv")$n_tests_regionwide); stopifnot(length(z) == 1); z }
N_PG_GROUPS   <- { z <- unique(rd("15_fourthcorner_pergroup.csv")$n_tests); stopifnot(length(z) == 1); z }
stopifnot(exists("HYP_ROWS"))
row_stats <- function(rows) {
  paste(vapply(rows, function(r) {
    parts <- c(if (!is.na(r$den)) sprintf("%s: %s", if (is.na(r$bio)) "" else "weighted by density", stat_str(r$den)),
               if (!is.na(r$bio)) sprintf("weighted by biomass: %s", stat_str(r$bio)))
    parts <- sub("^: ", "", parts)
    sprintf("%s (%s)", gsub("\n", " ", r$label), paste(parts, collapse = "; "))
  }, character(1)), collapse = ". ")
}
HYP_RULE <- "A row is drawn for every metric that met the screening criterion in at least one region under at least one weighting, and a density-biomass pair of the same metric is drawn whole, so the weighting that did not respond is seen beside the one that did (the two totals are such a pair, in two rows because their units differ); Rao's quadratic entropy is omitted because it duplicates functional dispersion (Table S10)."
HYP_ENCODING <- "Filled points and solid lines are weighted by density, open points and dashed lines weighted by biomass (for pollution tolerance, relative abundance and biomass); a trend (a generalized additive model smooth, mgcv with k = 3 under REML, with its 95% band) is drawn only where that metric and weighting met the criterion, and it describes the shape of the response without claiming linearity, since the test is a rank correlation. Color encodes physiographic region throughout the figure set; condition-class shapes are not drawn here because shape encodes weighting. Each row has its own y scale."
cap3 <- paste(
  "**Figure 3.** The metrics of hypotheses 1 and 2 that responded to the condition gradient, within region.",
  CRITERION, HYP_RULE, HYP_ENCODING, THRESH_NOTE, VR_TWO_BIOREGIONS,
  paste0("Rows, top to bottom, with the regions in which each met the criterion: ", row_stats(HYP_ROWS$fig3), "."),
  "Coefficients are in Tables 3a and 3b and every test in Figure 5.")
cap4 <- paste(
  "**Figure 4.** The metrics of hypothesis 3 that responded to the condition gradient, within region: overlap redundancy and the inequality measures under each weighting, the two totals, and the two density-biomass divergences.",
  CRITERION, HYP_RULE, HYP_ENCODING, THRESH_NOTE,
  paste0("Rows, top to bottom, with the regions in which each met the criterion: ", row_stats(HYP_ROWS$fig4), "."),
  "Coefficients are in Table 3b and every test in Figure 5.")
cap5 <- paste(
  "**Figure 5.** Every metric-by-weighting test against the condition gradient, within region, except Rao's quadratic entropy, which duplicates functional dispersion (Table S10).",
  CRITERION,
  "(a) Structure and totals; (b) functional diversity, redundancy, inequality and the density-biomass divergences: each point is Spearman's rho with its Fisher-z interval, filled where the test met the criterion, open where it did not, with the criterion drawn as the dashed lines.",
  "The intervals give the precision of each estimate and bound an individual coefficient, not the null; the detection bound is the power calculation (Table S14).",
  "The density-against-biomass difference of each trait state is Figure S2.",
  BOUNDED_NOTE)
capS2 <- paste(
  "**Figure S2.** The density-against-biomass difference of each trait state's fourth-corner statistic, by region.",
  sprintf("Each diamond is one trait state's fourth-corner statistic weighted by biomass minus its statistic weighted by density, for the %d states counted in each region (one per two-state group), and each gray bar is the 95%% percentile interval from %s resamples in which streams and taxa were resampled independently, with replacement, on each draw, drawn as its range over %d seeds; a filled diamond would mark a state whose difference met FDR <= %s within the region under all %d seeds, and an open diamond is a state not shown to differ between the weightings, which does not mean the two statistics are the same.",
          WD_PER_REGION, format(WD_RESAMPLES, big.mark = ","), N_SEEDS_FC, FC_FDR, N_SEEDS_FC),
  "A state that passes under one weighting and not the other in Figure 6 is not a tested difference between the weightings; this figure is that test, and no state met it.")

capS3 <- local({
  met <- rd("20_chironomid_functional_metrics.csv"); cnt <- rd("20_chironomid_fourthcorner_counts.csv"); sh <- rd("20_chironomid_share.csv")
  d <- dplyr::filter(met, Weighting == "Density", Metric %in% c("FDis_den", "RaoQ_den", "FR_RD_den"))
  rao <- dplyr::filter(met, Metric == "RaoQ_den")
  stopifnot("Figure S3 caption: a random set of four moved Rao's Q as far as removing the chironomid taxa; reword" = all(rao$random4_as_extreme == 0),
            "Figure S3 caption: a metric meeting the criterion without the four is not reversed in sign; reword" =
              all(with(dplyr::filter(d, met_removed), sign(rho_removed) != sign(rho_all))))
  paste(
    "**Figure S3.** Functional metrics and trait states weighted by density, with every taxon and without the four chironomid taxa (the tribes Chironomini and Tanytarsini and the subfamilies Orthocladiinae and Tanypodinae).",
    sprintf("The four made up %.0f%% to %.0f%% of the trait-set density (Table S26). (a) Each line joins a region's Spearman correlation between the metric, weighted by density, and the B-IBI with every taxon (black) to the same correlation without the four (gray); a filled point meets the screening criterion (|rho| of at least 0.30 and FDR of at most %.2f within its metric family and region, as in Table S10) and an open point does not. With every taxon, %d of the 9 tests met the criterion; without the four, %d did, in the opposite direction (Table S27).",
            min(sh[[3]][sh$Weighting == "Density"]), max(sh[[3]][sh$Weighting == "Density"]), FDR_MAX, sum(d$met_all), sum(d$met_removed)),
    sprintf("(b) The fourth-corner statistic of every trait state that passed under all five seeds (FDR of at most %.2f within its region-by-trait-group block, permutation model 6, %s permutations) in either analysis, one state per two-state group, with every taxon (black) and without the four (gray); points are offset vertically so neither hides the other. Weighted by density, %d of the %d states that passed with every taxon still passed without the four and %d passed only without them; weighted by biomass, %d of %d, and %d (Table S28).",
            FDR_MAX, format(N_PERM_FC, big.mark = ","), sum(cnt$passing_all[cnt$Weighting == "Density"]) - sum(cnt$lost[cnt$Weighting == "Density"]), sum(cnt$passing_all[cnt$Weighting == "Density"]), sum(cnt$gained[cnt$Weighting == "Density"]),
            sum(cnt$passing_all[cnt$Weighting == "Biomass"]) - sum(cnt$lost[cnt$Weighting == "Biomass"]), sum(cnt$passing_all[cnt$Weighting == "Biomass"]), sum(cnt$gained[cnt$Weighting == "Biomass"])),
    sprintf("Removing the four most abundant other taxa left Rao's quadratic entropy weighted by density associated with condition in every region (rho %s to %s), and no random set of four other taxa moved it as far as removing the four chironomid taxa did (Table S27).",
            sprintf("%+.3f", min(rao$rho_top4_other[rao$Weighting == "Density"])), sprintf("%+.3f", max(rao$rho_top4_other[rao$Weighting == "Density"]))))
})

# The two result clauses in the legacy Figure 6 caption are asserted against the tables of this run (27 Sep 2026):
# the per-group test on the fuzzy correspondence analysis admits no group, and no ring or box appears under biomass.
stopifnot("Figure 6 caption: the fuzzy-CA per-group test now admits a group; reword the caption" =
            all(dplyr::bind_rows(rd("15_fourthcorner_pergroup.csv"), rd("15_fourthcorner_pergroup_biomass.csv")) %>%
                  dplyr::filter(grepl("^fuzzy correspondence", q_ordination)) %>% dplyr::pull(n_seeds_pass) < N_SEEDS_FC),
          "Figure 6 caption: a ring or box now appears under biomass; reword the caption" =
            all(rd("15_fourthcorner_regionwide_biomass.csv")$n_seeds_pass_regionwide < N_SEEDS_FC) &&
            all(rd("15_fourthcorner_pergroup_biomass.csv") %>% dplyr::filter(grepl("^states standardized", q_ordination)) %>% dplyr::pull(n_seeds_pass) < N_SEEDS_FC))
cap6 <- paste(
  "**Figure 6.** Trait states associated with biological condition within their trait groups, and the stricter corrections each also met, by region, weighted by density (circles) and weighted by biomass (triangles).",
  sprintf("A filled mark is one trait state whose fourth-corner association with the index met FDR <= %s within its trait group under all %d seeds (%s permutations), an open mark a borderline state that met it under some seeds only, a ring a state that also met it across the region's %d states, and a box a trait group that met it in one test per group across the region's %d groups with states standardized; the same per-group test on the fuzzy correspondence analysis, the rule fixed first, admits no group.",
          FC_FDR, N_SEEDS_FC, format(N_PERM_FC, big.mark = ","), N_RW_STATES, N_PG_GROUPS),
  "Rings and boxes appear only in the columns weighted by density, and a state that passes under one weighting and not the other is not a tested difference between the weightings, which Figure S2 tests.",
  "One mark is drawn per two-state group, whose states are mirrors; statistics, seed ranges and carrier taxa are in Table S16.")

# Caption of the binned trait-space figure (capS1; not written, because the figure is not carried). Every number is
# read from the table that holds it: the trait set from step 11, each region's pool from step 16, r from step 11's
# dimensionality table.
S1_NTRAIT <- rd("11_community_vs_trait_set.csv")$n_trait_set
S1_POOL   <- { z <- rd("16_cross_region_pool_beta.csv")
               z <- dplyr::bind_rows(dplyr::transmute(z, Region = Region_A, taxa = taxa_A), dplyr::transmute(z, Region = Region_B, taxa = taxa_B)) %>% dplyr::distinct()
               stopifnot(!any(duplicated(z$Region))); setNames(z$taxa, z$Region) }
S1_PCT    <- function(rg) sprintf("%.1f", 100 * S1_POOL[[rg]] / S1_NTRAIT)
S1_R2     <- { z <- rd("11_traitspace_dimensionality.csv") %>% dplyr::filter(n_axes == 2); stopifnot(nrow(z) == 1); sprintf("%.2f", z$pearson_r) }
capS1 <- paste("**Figure S1.** Trait space occupied by sites in the bottom and top third of the within-region B-IBI gradient, weighted by density and by biomass, in each region.",
  sprintf("The backdrop separates the potential from the realized species pool: open circles are taxa present in the %d-taxon metacommunity but never collected in that region, filled circles are taxa collected there.", S1_NTRAIT),
  sprintf("The taxa collected in Valley and Ridge make up %s%% of the metacommunity, those in the Piedmont %s%% and those in the Coastal Plain %s%%, so part of the occupied trait space is unavailable to each region before condition is considered at all.",
          S1_PCT("Valley & Ridge"), S1_PCT("Piedmont"), S1_PCT("Coastal Plain")),
  "Colored points are site scores, hulls enclose the sites in each condition group, and the arrow runs from the low-condition to the high-condition centroid.",
  sprintf("This ordination is a low-dimensional aid, not evidence: two principal coordinate axes reproduce the underlying trait dissimilarities at r = %s (Table S6; Maire et al., 2015), so the inference rests on the functional metrics, which are computed from the full dissimilarity matrix.", S1_R2))

# capS2 above is the caption of the trait-state difference figure (figS4_traitstate_difference).

# Graphical abstract caption (6 Oct 2026). The three region coefficients it quotes are the divergence row of the
# screen; the caption claims the divergence met the criterion in every region, so the run halts if that stops being true.
ga_rho <- function(metric, region) {
  r <- screen %>% filter(Metric == metric, Region == region); stopifnot(nrow(r) == 1); sprintf("%+.2f", r$rho) }
stopifnot("Graphical abstract caption: the divergence no longer meets the criterion in every region; reword it" =
            all((screen %>% filter(Metric == "FR_DB_RD"))$Tier == "strong"),
          "Graphical abstract caption: the side that differs has changed (Valley & Ridge and Coastal Plain density, Piedmont biomass); reword it" =
            all((screen %>% filter(Metric == "FR_RD_den", Region %in% c("Valley & Ridge", "Coastal Plain")))$Tier == "strong") &&
            identical((screen %>% filter(Metric == "FR_RD_bio", Region == "Piedmont"))$Tier, "strong"))
# 6 Oct 2026: "the numerous taxa" -> "the high-density taxa". Caption text only.
capGA <- paste(
  "**Graphical abstract.** Trait sharing (functional redundancy) of stream macroinvertebrates weighted by density (individuals m-2) and by biomass (mg m-2) in the poorer and the better third of 30 streams in each of three Chesapeake Bay regions, ranked by biological condition (the regional benthic index of biotic integrity, B-IBI); the middle third is not drawn.",
  "Large points are medians of ten streams, small points single streams (density left, biomass right of each median), and the tinted area is the gap between the two weightings.",
  sprintf("In every region, biomass-weighted minus density-weighted sharing was higher where condition was better (Spearman rho = %s, %s, %s in Valley and Ridge, the Piedmont and the Coastal Plain; false discovery rate <= %s): where condition was poorer, the high-density taxa, which contribute most of the density, shared more traits than the large-bodied taxa, which contribute most of the biomass (proxy: body size is inferred from the weighting, not measured per taxon here).",
          ga_rho("FR_DB_RD", "Valley & Ridge"), ga_rho("FR_DB_RD", "Piedmont"), ga_rho("FR_DB_RD", "Coastal Plain"), format(FDR_MAX, nsmall = 2)),
  sprintf("The side that differed varied by region: density-weighted sharing in Valley and Ridge and the Coastal Plain (rho = %s, %s), biomass-weighted sharing in the Piedmont (rho = %s).",
          ga_rho("FR_RD_den", "Valley & Ridge"), ga_rho("FR_RD_den", "Coastal Plain"), ga_rho("FR_RD_bio", "Piedmont")))

# CAPTIONS. The manuscript carries the AUTHORED captions and is the authority for
#   their wording; it is not read here. This step generates captions whose every
#   coefficient is read live from the tables of this run, and CAPTIONS.md prints
#   them as the numbers check against the authored text.
GENERATED <- list(Figure_3 = cap3, Figure_4 = cap4, Figure_5 = cap5, Figure_6 = cap6, Figure_S2 = capS2, Figure_S3 = capS3, Graphical_abstract = capGA)   # legacy numbers; cap2 and capS1 are not written

emit_fig <- function(nm) {
  generated <- GENERATED[[nm]]
  if (is.null(generated))
    return(c(sprintf("**%s.** Authored caption in the manuscript; no coefficient in it comes from the screen, so none is generated here.",
                     sub("_", " ", nm)), ""))
  c("> **Pipeline caption, every coefficient read live from the tables of this run.**",
    "> Use it to check the authored caption in the manuscript; do not paste it into it.",
    ">", paste(">", generated), "")
}

cap <- c("# Chapter 2 -- figure and table captions", "",
  sprintf("Generated by pipeline step 18 on %s.", format(Sys.Date())),
  "",
  "**Figure captions are authored in the manuscript**, which is the single authority for",
  "their wording. Where the pipeline can produce the same coefficients, its caption is",
  "printed here as a numbers check, so a caption that has drifted from its figure is",
  "visible by reading the two side by side. Table captions below are generated from the asset map.",
  "",
  sprintf("Figures are drawn at %.2f in (174 mm) wide, %d dpi, %s, %.1f pt body type,",
          FIG_W_DOUBLE, FIG_DPI, FIG_FAMILY, FIG_BASE_PT),
  "with panel labels (a) (b) (c). Statistics are deliberately absent from the panels.",
  "",
  "Figure 1 is drawn by hand and has no pipeline source; Tables 1, 2, S2 and S3 are",
  "authored summaries and so appear in the manuscript but not in the table list below.",
  "", "---", "", "## Main figures", "")
for (nm in c("Figure_1", "Figure_2", "Figure_3", "Figure_4", "Figure_5", "Figure_6")) cap <- c(cap, emit_fig(nm))
cap <- c(cap, "## Graphical abstract (not numbered; J-B, 6 October 2026)", "", emit_fig("Graphical_abstract"))
cap <- c(cap, "## Main tables", "")
for (r in which(manifest$Placement == "main" & manifest$Kind == "table"))
  cap <- c(cap, sprintf("**%s.** %s. Carries: %s.", sub("_", " ", manifest$Asset[r]),
                        manifest$Title[r], manifest$Carries[r]), "")
cap <- c(cap, "## Supplementary figures", "")
for (nm in c("Figure_S2", "Figure_S3")) cap <- c(cap, emit_fig(nm))
cap <- c(cap, "## Supplementary tables", "")
for (r in which(manifest$Placement == "supplementary" & manifest$Kind == "table"))
  cap <- c(cap, sprintf("**%s.** %s. Carries: %s.", sub("_", " ", manifest$Asset[r]),
                        manifest$Title[r], manifest$Carries[r]), "")
cap <- c(cap, "## Record tables (the data behind a figure; not printed)", "")
for (r in which(manifest$Placement == "record" & manifest$Kind == "table"))
  cap <- c(cap, sprintf("**%s.** %s. Carries: %s.", sub("_", " ", manifest$Asset[r]),
                        manifest$Title[r], manifest$Carries[r]), "")
writeLines(cap, file.path(MS, "CAPTIONS.md"))

# WHICH SOURCE EACH CAPTION CAME FROM. A reader of the manifest
# should not have to open CAPTIONS.md to learn whether a caption is authored
# prose or pipeline output, because the two are edited in different places and
# only one of them is safe to edit by hand.
cap_source <- function(a) {
  if (!identical(manifest$Kind[manifest$Asset == a][1], "figure"))
    return("generated -- asset map; the manuscript carries the authored version")
  if (!is.null(GENERATED[[a]])) "authored in the manuscript; pipeline caption with live coefficients in CAPTIONS.md"
  else "authored in the manuscript; no pipeline caption"
}
manifest$Caption <- vapply(manifest$Asset, cap_source, character(1))
manifest$Source  <- ifelse(manifest$Asset == "Figure_1", "hand-drawn",
                    ifelse(manifest$Asset == "Graphical_abstract", "pipeline tables; Manuscripts/Graphical_abstract_2026-10-06/scripts/52_J_B_slope_pairs.py", "pipeline"))
write_csv(manifest, file.path(MS, "ASSET_MANIFEST.csv"))
say("  captions: %d pipeline captions written to CAPTIONS.md; the authored captions are in the manuscript", length(GENERATED))

say("  main assets       : %d figures, %d tables",
    sum(manifest$Placement == "main" & manifest$Kind == "figure"),
    sum(manifest$Placement == "main" & manifest$Kind == "table"))
say("  supplementary     : %d figures, %d tables",
    sum(manifest$Placement == "supplementary" & manifest$Kind == "figure"),
    sum(manifest$Placement == "supplementary" & manifest$Kind == "table"))
say("  written to output/manuscript/ with ASSET_MANIFEST.csv and CAPTIONS.md")
say("STEP 18 complete")
