# Data dictionary

Every file under `data/`, what it is, where its values come from, which script reads or writes it, and every column. Roles: **input** is a shipped record the pipeline reads; **derived** is written by `R/00c`, `R/00a` or `R/00b` from the inputs and shipped so that `run_all.R` can run alone; **documentation** is read by no step. Paths are relative to `data/`, which mirrors the project's own data folder, so the scripts run unchanged.

Units throughout: density in individuals per square meter, biomass in milligrams dry mass per square meter, the Chessie B-IBI on its 0 to 100 scale. Sites are the 90 wadeable streams of the study; regions are Valley and Ridge, Piedmont and Coastal Plain.

## Inputs

### `CONUS_Ancillary_Taxonomy.csv`

Order, family and genus of every genus in the CONUS trait database, with ITIS accepted names.

- Role: input
- Read by: R/00c_resolve_parent_rows.R (lineage of genus rows)
- Source: Twardochleb et al. (2021), CONUS trait database ancillary taxonomy
- One row is: one genus
- Rows: 5373; columns: 7

| Column | Description |
|---|---|
| `Submitted_name` | Name as submitted |
| `Accepted_name` | ITIS accepted name |
| `Accepted_TSN` | ITIS serial number |
| `Order` | Order |
| `Family` | Family |
| `Genus` | Genus |
| `Species` | Species where the record is a species, else NA |

### `H1_HBI_taxon_contributions.csv`

Every site by taxon with its family tolerance value and its weights and contributions to the site's HBI under the three weightings.

- Role: input
- Read by: R/00_config.R (IN$hbi_tv): R/09b_taxon_tolerance.R
- Source: Upstream Hilsenhoff Biotic Index computation with the Chessie B-IBI family tolerance values
- One row is: one site by taxon
- Rows: 4041; columns: 14

| Column | Description |
|---|---|
| `Typology` | Study region coded as the project's typology number: 1 = Valley and Ridge (Shenandoah), 2 = Coastal Plain (Delmarva), 3 = Piedmont (southeastern Pennsylvania and Maryland) |
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `TAXANAME` | Taxon name at the finest rank identified (genus where possible, otherwise family, subfamily, tribe or a higher rank); the key that joins the community table, the matrices and the trait matrix |
| `Family` | Family, upper case |
| `TV` | Family tolerance value, 0 to 10, NA where none |
| `abund_w` | Total corrected abundance |
| `den_w` | Total density |
| `bio_w` | Total biomass |
| `rel_abund` | Share of the site's abundance |
| `rel_density` | Share of the site's density |
| `rel_biomass` | Share of the site's biomass |
| `HBI_rel_abund_contribution` | TV times rel_abund |
| `HBI_density_contribution` | TV times rel_density |
| `HBI_biomass_contribution` | TV times rel_biomass |

### `Macros_Typ-1_ChessieB-IBI_Scores_Metrics.csv`

Output of the Chessie B-IBI scoring workflow of the Interstate Commission on the Potomac River Basin (Smith et al. 2017) run on the Valley and Ridge collections: the observed score, and the median, mean and rating over 100 rarefied subsamples of 100 individuals. The pipeline uses median_score as the biological condition gradient and median_rating as the class; step 01 stops unless every site's value in H1_BIBI_site_with_classes_MEDIAN.csv equals its median_score here. The remaining columns are the index's component metrics as the workflow names them.

- Role: input
- Read by: R/01_load_and_verify.R (provenance gate); the project's class-file builder, build_bibi_site_classes.R (upstream, not shipped)
- Source: ICPRB Chessie B-IBI scoring workflow on the project's collections
- One row is: one site
- Rows: 30; columns: 75

| Column | Description |
|---|---|
| `unique_id` | ICPRB sample identifier |
| `icprb_bioregion_id` | ICPRB bioregion used to score the site (sgv, srv, unp, mac) |
| `site_ID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `site_info` | Site name |
| `COMID` | NHDPlus reach identifier |
| `date` | Sample date |
| `chessie_bibi` | Observed Chessie B-IBI score, 0 to 100, on the full sample |
| `rating` | Narrative rating of the observed score |
| `median_score` | Median score over 100 rarefactions to 100 individuals; the gradient the paper uses |
| `median_rating` | Narrative rating of the median score; the class the paper uses |
| `median_percent` | Percent of rarefactions giving the median rating |
| `mean_score` | Mean score over the rarefactions |
| `st_dev` | Standard deviation of the score over the rarefactions |
| `mean_rating` | Rating of the mean score |
| `mean_percent` | Percent of rarefactions giving the mean rating |
| `aspt_mod`, `pct_gastro_oligo`, `pct_diptera`, `gold`, `hbi`, `pct_arthropoda` and 54 more (pattern `^(pct_|rich_|aspt_mod|gold|hbi|margalefs|pielou)`) | Component metric of the Chessie B-IBI as the ICPRB workflow names it (Smith et al. 2017, Table 5) |

### `Macros_Typ-2_ChessieB-IBI_Scores_Metrics.csv`

Output of the Chessie B-IBI scoring workflow of the Interstate Commission on the Potomac River Basin (Smith et al. 2017) run on the Coastal Plain collections: the observed score, and the median, mean and rating over 100 rarefied subsamples of 100 individuals. The pipeline uses median_score as the biological condition gradient and median_rating as the class; step 01 stops unless every site's value in H1_BIBI_site_with_classes_MEDIAN.csv equals its median_score here. The remaining columns are the index's component metrics as the workflow names them.

- Role: input
- Read by: R/01_load_and_verify.R (provenance gate); the project's class-file builder, build_bibi_site_classes.R (upstream, not shipped)
- Source: ICPRB Chessie B-IBI scoring workflow on the project's collections
- One row is: one site
- Rows: 30; columns: 75

| Column | Description |
|---|---|
| `unique_id` | ICPRB sample identifier |
| `icprb_bioregion_id` | ICPRB bioregion used to score the site (sgv, srv, unp, mac) |
| `site_id` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `station_id` | Station identifier |
| `HUC12` | Twelve-digit hydrologic unit |
| `sample_date` | Sample date |
| `chessie_bibi` | Observed Chessie B-IBI score, 0 to 100, on the full sample |
| `rating` | Narrative rating of the observed score |
| `median_score` | Median score over 100 rarefactions to 100 individuals; the gradient the paper uses |
| `median_rating` | Narrative rating of the median score; the class the paper uses |
| `median_percent` | Percent of rarefactions giving the median rating |
| `mean_score` | Mean score over the rarefactions |
| `st_dev` | Standard deviation of the score over the rarefactions |
| `mean_rating` | Rating of the mean score |
| `mean_percent` | Percent of rarefactions giving the mean rating |
| `aspt_mod`, `pct_gastro_oligo`, `pct_diptera`, `gold`, `hbi`, `pct_arthropoda` and 54 more (pattern `^(pct_|rich_|aspt_mod|gold|hbi|margalefs|pielou)`) | Component metric of the Chessie B-IBI as the ICPRB workflow names it (Smith et al. 2017, Table 5) |

### `Regional Analysis/env/tables/H2_gradient_scores_by_typology.csv`

The land-use gradient along which the sites were selected: the first principal component of a within-region land-use PCA, one row per site. It labels the study design in Table S1 and enters no test.

- Role: input
- Read by: R/18_manuscript_assets.R (Table S1 part A's site-selection gradient column only; no analysis uses it)
- Source: Within-region principal component analysis of catchment land use, upstream
- One row is: one site
- Rows: 90; columns: 8

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Typology` | Study region coded as the project's typology number: 1 = Valley and Ridge (Shenandoah), 2 = Coastal Plain (Delmarva), 3 = Piedmont (southeastern Pennsylvania and Maryland) |
| `PC1` | Site score on the first land-use axis |
| `PC2` | Second axis |
| `PC3` | Third axis |
| `PC1_flipped` | PC1 with its sign set so that higher means more agricultural or developed |
| `PC1_sig_cor` | Whether PC1 correlated with the index in the upstream check |
| `ag_grad_lu` | The selection gradient as used: PC1_flipped |

### `Regional Analysis/tables/H1/H1_BIBI_site_with_classes_MEDIAN.csv`

The biological condition gradient: each site's median Chessie B-IBI score and its ICPRB rating, with the bioregion whose cut-points scored it.

- Role: input
- Read by: R/00_config.R (IN$bibi under CH2_BIBI_MEDIAN=1): steps 01, 19 and every step through the site table
- Source: Written from the three ICPRB files by Data/Data/build_bibi_site_classes.R (not shipped; it also reads the ICPRB scoring scripts for the cut-points); never edited by hand
- One row is: one site
- Rows: 90; columns: 6

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Typology` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `BIBI` | Median Chessie B-IBI over 100 rarefactions, 0 to 100, equal to the ICPRB median_score |
| `BIBI_Class` | ICPRB median_rating: Very poor, Poor, Fair, Good or Excellent |
| `Region` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `Bioregion` | ICPRB bioregion whose cut-points scored the site: sgv or srv (Valley and Ridge), unp (Piedmont), mac (Coastal Plain) |

### `Regional Analysis/tables/H1/H1_BIBI_thresholds_used.csv`

The four class boundaries of each bioregion, and how many of the region's sites each bioregion scored.

- Role: input
- Read by: R/18_manuscript_assets.R (the rating-class cut-point columns of Table S1 part A)
- Source: ICPRB scoring scripts, section 5 (Smith et al. 2017, Table 5), through build_bibi_site_classes.R
- One row is: one region by boundary by bioregion
- Rows: 16; columns: 5

| Column | Description |
|---|---|
| `Typology` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `Label` | Boundary: VeryPoor–Poor, Poor–Fair, Fair–Good or Good–Excellent |
| `Score` | Score at the boundary, 0 to 100 |
| `Bioregion` | ICPRB bioregion |
| `n_sites` | Sites of the region scored under this bioregion |

### `Regional Analysis/tables/H1/H1_site_summary_macro_metrics_twardochleb_unresolved_2026-10-06.csv`

The unresolved site table after R/00a rebuilt its trait columns from the Twardochleb fuzzy matrix and dropped the legacy columns.

- Role: input
- Read by: R/00c_resolve_parent_rows.R (carries the community-weighted means for R/00a)
- Source: R/00a run on the unresolved data, 22 September 2026
- One row is: one site
- Rows: 90; columns: 118

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Typology` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `FRic_den_fuzPool` | Functional richness computed in the trait space of the whole taxon pool (all three regions) under the fuzzy coding; read by R/00b for space-filling redundancy |
| `FRic_den` | Functional richness: convex hull volume of the site's trait set in the ten-axis PCoA of the fuzzy group-weighted Gower distance (FD::dbFD, Cailliez correction); presence-based, one value per site |
| `FEve_den`, `FDiv_den`, `FDis_den`, `RaoQ_den`, `FEve_bio`, `FDiv_bio` and 2 more (pattern `^(FEve|FDiv|FDis|RaoQ)_(den|bio)$`) | Functional evenness, divergence, dispersion or Rao's quadratic entropy (FD::dbFD on the fuzzy proxy trait matrix with block-balanced Gower weights, Cailliez correction, ten axes), weighted by mean density (den) or mean biomass (bio) over the processed replicates; Villeger et al. 2008, Laliberte and Legendre 2010, Botta-Dukat 2005 |
| `Simpson_den` | Simpson diversity 1 - sum(p_i^2) on mean density, over the taxa with a trait record |
| `Simpson_bio` | Simpson diversity 1 - sum(p_i^2) on mean biomass, over the taxa with a trait record |
| `Richness_FD` | Number of taxa with a trait record and a nonzero total at the site (the trait set) |
| `RaoQ_gower_den` | Rao's quadratic entropy on the binary unweighted Gower distance, density weights; the geometry the upstream table used for redundancy, replaced by R/00b's fuzzy geometry |
| `RaoQ_gower_bio` | As RaoQ_gower_den, biomass weights |
| `FR_RD_den` | Overlap redundancy weighted by density, 1 - Q/D, with Q Rao's quadratic entropy on the fuzzy group-weighted Gower distance and D Simpson diversity (de Bello et al. 2007; Ricotta et al. 2016) |
| `FR_RD_bio` | Overlap redundancy weighted by biomass, 1 - Q/D |
| `FR_SF_den` | Space-filling redundancy: Richness_FD divided by FRic_den_fuzPool (Mouillot et al. 2014) |
| `FR_SF_den_log10` | log10 of FR_SF_den; the screened form |
| `FR_DB_RD` | Redundancy divergence: FR_RD_bio minus FR_RD_den |
| `Richness_alpha` | Taxonomic richness: number of taxa present at the site, at the finest rank identified |
| `Shannon` | Shannon diversity H' on mean density over the processed replicates, natural logarithm, over every taxon |
| `Total_Density` | Total density, individuals per square meter: the median over the site's processed replicate samples of each replicate's total (R/00d, from Step2_SiteReplicateMedians.csv; the lab method M3). On the unresolved tables, before R/00d, the mean over the processed replicates |
| `Total_Biomass` | Total biomass, milligrams dry mass per square meter: the median over the processed replicates of each replicate's total (R/00d; M3). On the unresolved tables, the mean over the processed replicates |
| `BIBI_score` | Chessie B-IBI score of the site as the upstream site table carried it; the pipeline reads the index from H1_BIBI_site_with_classes_MEDIAN.csv instead, so this column is carried, never used |
| `BIBI_Class` | Narrative rating carried from the upstream site table; not used by the pipeline (see BIBI_score) |
| `HBI_taxa_with_TV` | Number of taxa at the site with a tolerance value |
| `HBI_taxa_total` | Number of taxa at the site |
| `HBI_rel_abund` | Hilsenhoff Biotic Index weighted by relative abundance (the conventional form; family tolerance values from the Chessie B-IBI, Buchanan et al. 2011); the median over the processed replicates of each replicate's index after R/00d, the index on the mean taxon values before it |
| `HBI_den` | Hilsenhoff Biotic Index weighted by density; identical to HBI_rel_abund within a site and not screened |
| `HBI_bio` | Hilsenhoff Biotic Index weighted by biomass (median over the processed replicates after R/00d) |
| `HBI_coverage_pct` | Percent of abundance in taxa with a tolerance value |
| `HBI_density_coverage_pct` | Percent of density in taxa with a tolerance value |
| `HBI_biomass_coverage_pct` | Percent of biomass in taxa with a tolerance value |
| `Flying_Strong_CWM_den`, `Flying_Weak_CWM_den`, `Flying_Strong_CWM_bio`, `Flying_Weak_CWM_bio` | Adult flying strength: strong or weak. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Emerge_Fall_CWM_den`, `Emerge_Spring_CWM_den`, `Emerge_Summer_CWM_den`, `Emerge_Winter_CWM_den`, `Emerge_Fall_CWM_bio`, `Emerge_Spring_CWM_bio` and 2 more (pattern `^Emerge_(Fall|Spring|Summer|Winter)_CWM_(den|bio)$`) | Emergence season: fall, spring, summer or winter. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Synchro_Poorly_CWM_den`, `Synchro_Well_CWM_den`, `Synchro_Poorly_CWM_bio`, `Synchro_Well_CWM_bio` | Synchronization of emergence: poorly or well synchronized. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Feed_CF_CWM_den`, `Feed_CG_CWM_den`, `Feed_HB_CWM_den`, `Feed_PA_CWM_den`, `Feed_PR_CWM_den`, `Feed_SH_CWM_den` and 6 more (pattern `^Feed_(CF|CG|HB|PA|PR|SH)_CWM_(den|bio)$`) | Feeding style: CF collector-filterer, CG collector-gatherer, HB herbivore, PA parasite, PR predator, SH shredder. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Disperse_High_CWM_den`, `Disperse_Low_CWM_den`, `Disperse_High_CWM_bio`, `Disperse_Low_CWM_bio` | Female dispersal: high or low. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Habit_Burrower_CWM_den`, `Habit_Climber_CWM_den`, `Habit_Clinger_CWM_den`, `Habit_Crawler_CWM_den`, `Habit_Planktonic_CWM_den`, `Habit_Skater_CWM_den` and 10 more (pattern `^Habit_(Burrower|Climber|Clinger|Crawler|Planktonic|Skater|Sprawler|Swimmer)_CWM_(den|bio)$`) | Habit. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Size_Large_CWM_den`, `Size_Medium_CWM_den`, `Size_Small_CWM_den`, `Size_Large_CWM_bio`, `Size_Medium_CWM_bio`, `Size_Small_CWM_bio` | Maximum body size class. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Respiration_Gills_CWM_den`, `Respiration_Plastron_CWM_den`, `Respiration_Tegument_CWM_den`, `Respiration_Gills_CWM_bio`, `Respiration_Plastron_CWM_bio`, `Respiration_Tegument_CWM_bio` | Respiration mode; Plastron covers plastron and spiracle (Twardochleb's Plastron_spiracle) and is never merged with gills. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Rheophily_Both_CWM_den`, `Rheophily_Depositional_CWM_den`, `Rheophily_Erosional_CWM_den`, `Rheophily_Both_CWM_bio`, `Rheophily_Depositional_CWM_bio`, `Rheophily_Erosional_CWM_bio` | Rheophily: depositional, erosional or both. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Thermal_ColdEury_CWM_den`, `Thermal_ColdSteno_CWM_den`, `Thermal_CoolEury_CWM_den`, `Thermal_HotEury_CWM_den`, `Thermal_WarmEury_CWM_den`, `Thermal_ColdEury_CWM_bio` and 4 more (pattern `^Thermal_(ColdEury|ColdSteno|CoolEury|HotEury|WarmEury)_CWM_(den|bio)$`) | Thermal preference: cold stenothermal, cold-cool eurythermal, cool-warm eurythermal, warm eurythermal, hot eurythermal. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Voltinism_Multivoltine_CWM_den`, `Voltinism_Semivoltine_CWM_den`, `Voltinism_Univoltine_CWM_den`, `Voltinism_Multivoltine_CWM_bio`, `Voltinism_Semivoltine_CWM_bio`, `Voltinism_Univoltine_CWM_bio` | Voltinism. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |

### `Regional Analysis/tables/H1/H1_site_summary_macro_metrics_unresolved_2026-10-06.csv`

The site table before parent rows were resolved: every structural, tolerance, functional diversity and redundancy metric per site as the upstream scripts computed them, plus legacy columns the pipeline no longer reads. R/00c recomputes every metric it later rewrites from this table's community data first, and stops unless each reproduces (the provenance gate), so the definitions applied to the resolved data are the ones that produced this table.

- Role: input
- Read by: R/00c_resolve_parent_rows.R (provenance gate, and the columns it carries)
- Source: Upstream site table (Step3), copy taken before R/00c first ran
- One row is: one site
- Rows: 90; columns: 163

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Typology` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `FRic_den_binPool`, `FEve_den_binPool`, `FDiv_den_binPool`, `FDis_den_binPool`, `RaoQ_den_binPool`, `FRic_bio_binPool` and 33 more (pattern `^(FRic|FEve|FDiv|FDis|RaoQ)_(den|bio)_(binPool|fuzPool|binReg|fuzReg)$`) | Legacy variant of the metric: binary (bin) or fuzzy (fuz) coding, in the pooled (Pool) or regional (Reg) trait space; dropped by R/00a except FRic_den_fuzPool |
| `FRic_den_fuzPool` | Functional richness computed in the trait space of the whole taxon pool (all three regions) under the fuzzy coding; read by R/00b for space-filling redundancy |
| `FRic_den` | Functional richness: convex hull volume of the site's trait set in the ten-axis PCoA of the fuzzy group-weighted Gower distance (FD::dbFD, Cailliez correction); presence-based, one value per site |
| `FEve_den`, `FDiv_den`, `FDis_den`, `RaoQ_den`, `FEve_bio`, `FDiv_bio` and 2 more (pattern `^(FEve|FDiv|FDis|RaoQ)_(den|bio)$`) | Functional evenness, divergence, dispersion or Rao's quadratic entropy (FD::dbFD on the fuzzy proxy trait matrix with block-balanced Gower weights, Cailliez correction, ten axes), weighted by mean density (den) or mean biomass (bio) over the processed replicates; Villeger et al. 2008, Laliberte and Legendre 2010, Botta-Dukat 2005 |
| `Simpson_den` | Simpson diversity 1 - sum(p_i^2) on mean density, over the taxa with a trait record |
| `Simpson_bio` | Simpson diversity 1 - sum(p_i^2) on mean biomass, over the taxa with a trait record |
| `Richness_FD` | Number of taxa with a trait record and a nonzero total at the site (the trait set) |
| `RaoQ_gower_den` | Rao's quadratic entropy on the binary unweighted Gower distance, density weights; the geometry the upstream table used for redundancy, replaced by R/00b's fuzzy geometry |
| `RaoQ_gower_bio` | As RaoQ_gower_den, biomass weights |
| `FR_RD_den` | Overlap redundancy weighted by density, 1 - Q/D, with Q Rao's quadratic entropy on the fuzzy group-weighted Gower distance and D Simpson diversity (de Bello et al. 2007; Ricotta et al. 2016) |
| `FR_RD_bio` | Overlap redundancy weighted by biomass, 1 - Q/D |
| `FR_SF_den` | Space-filling redundancy: Richness_FD divided by FRic_den_fuzPool (Mouillot et al. 2014) |
| `FR_SF_den_log10` | log10 of FR_SF_den; the screened form |
| `FR_DB_RD` | Redundancy divergence: FR_RD_bio minus FR_RD_den |
| `Richness_alpha` | Taxonomic richness: number of taxa present at the site, at the finest rank identified |
| `Shannon` | Shannon diversity H' on mean density over the processed replicates, natural logarithm, over every taxon |
| `Total_Density` | Total density, individuals per square meter: the median over the site's processed replicate samples of each replicate's total (R/00d, from Step2_SiteReplicateMedians.csv; the lab method M3). On the unresolved tables, before R/00d, the mean over the processed replicates |
| `Total_Biomass` | Total biomass, milligrams dry mass per square meter: the median over the processed replicates of each replicate's total (R/00d; M3). On the unresolved tables, the mean over the processed replicates |
| `BIBI_score` | Chessie B-IBI score of the site as the upstream site table carried it; the pipeline reads the index from H1_BIBI_site_with_classes_MEDIAN.csv instead, so this column is carried, never used |
| `BIBI_Class` | Narrative rating carried from the upstream site table; not used by the pipeline (see BIBI_score) |
| `Disperse_High.dispersal_CWM_den`, `Disperse_Low.dispersal_CWM_den`, `Emerge_Spring_CWM_den`, `Emerge_Summer_CWM_den`, `Emerge_Winter_CWM_den`, `Feed_CF_Collector-filterer_CWM_den` and 70 more (pattern `^(Flying|Emerge|Synchro|Feed|Disperse|Habit|Size|Respiration|Rheophily|Thermal|Voltinism)_.*_CWM_(den|bio)$`) | Legacy community-weighted mean from the binary trait file (a taxon with several states counted in full in each); dropped by R/00a |
| `CWM_Response_Density`, `CWM_Response_Biomass`, `CWM_Effect_Density`, `CWM_Effect_Biomass`, `CWM_Linking_Density`, `CWM_Linking_Biomass` and 6 more (pattern `^(CWM|REL)_(Effect|Linking|Response)_(Density|Biomass|den|bio)$`) | Legacy aggregate of the binary trait file by trait role (effect, linking, response); dropped by R/00a |
| `HBI_taxa_with_TV` | Number of taxa at the site with a tolerance value |
| `HBI_taxa_total` | Number of taxa at the site |
| `HBI_rel_abund` | Hilsenhoff Biotic Index weighted by relative abundance (the conventional form; family tolerance values from the Chessie B-IBI, Buchanan et al. 2011); the median over the processed replicates of each replicate's index after R/00d, the index on the mean taxon values before it |
| `HBI_den` | Hilsenhoff Biotic Index weighted by density; identical to HBI_rel_abund within a site and not screened |
| `HBI_bio` | Hilsenhoff Biotic Index weighted by biomass (median over the processed replicates after R/00d) |
| `HBI_coverage_pct` | Percent of abundance in taxa with a tolerance value |
| `HBI_density_coverage_pct` | Percent of density in taxa with a tolerance value |
| `HBI_biomass_coverage_pct` | Percent of biomass in taxa with a tolerance value |

### `SiteTaxa_Density_Biomass_Transformed_unresolved_2026-10-06.csv`

The community table before parent rows are resolved: one row per site by taxon. Density and biomass follow the lab's method (James et al. 2022; the default method M3+M0 of the data chain, 5 October 2026): each replicate sample's count, after the subsample expansion, is divided by its 0.25 square meter area, and a site-by-taxon value is the mean over the site's processed replicates, zeros included, with no weighting by the share of replicates in which the taxon occurred. Replicates that were lost are missing, not zero, so the divisor is the number of replicates processed (five at 86 sites, four at SPOU, KNIT and THOM, three at SOUT). Under this method the Mean_, Median_ and Total_ columns of a taxon hold the same value. The suffix _unresolved_2026-10-06 marks the copy taken before R/00c first ran; R/00c reads only this copy, so it can run again without resolving twice.

- Role: input
- Read by: R/00c_resolve_parent_rows.R
- Source: Field collections of the U.S. Geological Survey best management practice study, 2021 to 2023, as published in the three USGS data releases; processed upstream from raw counts and body lengths (length-mass regressions of Benke et al. 1999 and others) to per-sample density and biomass
- One row is: one site by taxon
- Rows: 4041; columns: 26

| Column | Description |
|---|---|
| `Typology` | Study region coded as the project's typology number: 1 = Valley and Ridge (Shenandoah), 2 = Coastal Plain (Delmarva), 3 = Piedmont (southeastern Pennsylvania and Maryland) |
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `TAXANAME` | Taxon name at the finest rank identified (genus where possible, otherwise family, subfamily, tribe or a higher rank); the key that joins the community table, the matrices and the trait matrix |
| `Family` | Family of the taxon, upper case, as recorded in the field data |
| `Mean_Corrected_Abundance` | Mean over the processed replicates of the count per sample, zeros included |
| `Median_Corrected_Abundance` | Same value as Mean_Corrected_Abundance under the lab method (the mean over the processed replicates) |
| `Total_Corrected_Abundance` | Same value under the lab method: the sum of the counts divided by the number of processed replicates |
| `Mean_Density_per_m2` | Mean over the processed replicates of density, individuals per square meter (count divided by 0.25 square meter), zeros included |
| `Median_Density_per_m2` | Same value as Mean_Density_per_m2 under the lab method; the density the pipeline uses |
| `Total_Density_per_m2` | Same value under the lab method: the sum of per-sample densities divided by the number of processed replicates |
| `Mean_Biomass_per_m2` | Mean over the processed replicates of biomass, milligrams dry mass per square meter, zeros included |
| `Median_Biomass_per_m2` | Same value as Mean_Biomass_per_m2 under the lab method; the biomass the pipeline uses |
| `Total_Biomass_per_m2` | Same value under the lab method: the sum of per-sample biomass divided by the number of processed replicates |
| `TAXANAME_raw` | Taxon name as it stood before name cleaning |
| `Mean_Corrected_Abundance_log`, `Median_Corrected_Abundance_log`, `Mean_Density_log`, `Mean_Biomass_log`, `Median_Density_log`, `Median_Biomass_log` | Natural logarithm of the same-named column plus one; carried, not used by the pipeline |
| `Mean_Corrected_Abundance_sqrt`, `Median_Corrected_Abundance_sqrt`, `Mean_Density_sqrt`, `Mean_Biomass_sqrt`, `Median_Density_sqrt`, `Median_Biomass_sqrt` | Square root of the same-named column; carried, not used by the pipeline |

### `Step2_SiteReplicateMedians.csv`

Site totals and the site Hilsenhoff Biotic Index by the lab's method (James et al. 2022): the median over the site's processed replicate samples of each replicate's value, with the mean beside it. R/00d writes the medians into Total_Density, Total_Biomass, HBI_rel_abund, HBI_den and HBI_bio of the two site tables and keeps the mean versions in the columns ending _mean, after checking that the mean totals equal the site tables' totals. Without this file R/00d changes nothing, and the site totals stay the mean version.

- Role: input
- Read by: R/00d_site_replicate_medians.R
- Source: Written by the data chain's Step 2 (Step2_Macroinvertebrate_Metric_Transformations, not shipped) under the default method M3+M0, from each replicate sample's density, biomass and Hilsenhoff index (6 October 2026)
- One row is: one site
- Rows: 90; columns: 11

| Column | Description |
|---|---|
| `Typology` | Study region coded as the project's typology number: 1 = Valley and Ridge (Shenandoah), 2 = Coastal Plain (Delmarva), 3 = Piedmont (southeastern Pennsylvania and Maryland) |
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `n_processed` | Replicate samples processed at the site (a lost replicate is missing, not zero) |
| `Total_Density_median` | Median over the processed replicates of the replicate's total density, individuals per square meter |
| `Total_Density_mean` | Mean over the processed replicates of the replicate's total density |
| `Total_Biomass_median` | Median over the processed replicates of the replicate's total biomass, milligrams dry mass per square meter |
| `Total_Biomass_mean` | Mean over the processed replicates of the replicate's total biomass |
| `HBI_rel_abund_median` | Median over the processed replicates of the replicate's Hilsenhoff Biotic Index weighted by relative abundance |
| `HBI_den_median` | Median of the replicate's index weighted by density (equal to the relative-abundance index) |
| `HBI_bio_median` | Median of the replicate's index weighted by biomass |
| `density_method` | Density method of the run that wrote the file; R/00d requires M3M0, the lab method |

### `Trait_Data/closest_kin_itis.csv`

Family, subfamily and tribe of every genus in the trait matrix, from ITIS.

- Role: input
- Read by: R/00c_resolve_parent_rows.R (to find the children of a subfamily or tribe row)
- Source: Integrated Taxonomic Information System (ITIS), getFullHierarchyFromTSN, accessed 22 September 2026
- One row is: one genus
- Rows: 198; columns: 7

| Column | Description |
|---|---|
| `Genus` | Genus |
| `Family` | Family (ITIS) |
| `Subfamily` | Subfamily (ITIS), blank where none |
| `Tribe` | Tribe (ITIS), blank where none |
| `Accepted_TSN` | ITIS taxonomic serial number of the accepted name |
| `source` | The ITIS call and access date |
| `status` | ok, or a note where ITIS returned nothing and the placement was set by hand |

### `Typ3_ChessieB-IBI_Scores_Metrics.csv`

Output of the Chessie B-IBI scoring workflow of the Interstate Commission on the Potomac River Basin (Smith et al. 2017) run on the Piedmont collections: the observed score, and the median, mean and rating over 100 rarefied subsamples of 100 individuals. The pipeline uses median_score as the biological condition gradient and median_rating as the class; step 01 stops unless every site's value in H1_BIBI_site_with_classes_MEDIAN.csv equals its median_score here. The remaining columns are the index's component metrics as the workflow names them.

- Role: input
- Read by: R/01_load_and_verify.R (provenance gate); the project's class-file builder, build_bibi_site_classes.R (upstream, not shipped)
- Source: ICPRB Chessie B-IBI scoring workflow on the project's collections
- One row is: one site
- Rows: 30; columns: 74

| Column | Description |
|---|---|
| `unique_id` | ICPRB sample identifier |
| `icprb_bioregion_id` | ICPRB bioregion used to score the site (sgv, srv, unp, mac) |
| `SITE_ID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `HUC12` | Twelve-digit hydrologic unit |
| `SAMPLE_DATE` | Sample date |
| `chessie_bibi` | Observed Chessie B-IBI score, 0 to 100, on the full sample |
| `rating100` | Narrative rating of the observed score |
| `median_score` | Median score over 100 rarefactions to 100 individuals; the gradient the paper uses |
| `median_rating` | Narrative rating of the median score; the class the paper uses |
| `median_percent` | Percent of rarefactions giving the median rating |
| `mean_score` | Mean score over the rarefactions |
| `st_dev` | Standard deviation of the score over the rarefactions |
| `mean_rating` | Rating of the mean score |
| `mean_percent` | Percent of rarefactions giving the mean rating |
| `aspt_mod`, `pct_gastro_oligo`, `pct_diptera`, `gold`, `hbi`, `pct_arthropoda` and 54 more (pattern `^(pct_|rich_|aspt_mod|gold|hbi|margalefs|pielou)`) | Component metric of the Chessie B-IBI as the ICPRB workflow names it (Smith et al. 2017, Table 5) |

### `bibi_tv_family_lookup.csv`

One tolerance value per family, lower case, with how many observations carried it and in which regions.

- Role: input
- Read by: R/00c_resolve_parent_rows.R (rebuilding the tolerance metrics)
- Source: Chessie B-IBI family tolerance values (Buchanan et al. 2011) as applied in the three ICPRB scoring runs
- One row is: one family
- Rows: 83; columns: 5

| Column | Description |
|---|---|
| `family` | Family, lower case |
| `bibi_tv` | Tolerance value, 0 to 10 |
| `n_obs` | Observations carrying it |
| `source_typologies` | Regions in which it occurred (T1, T2, T3) |
| `consistency_note` | consistent, or a note where the runs disagreed |

### `bio_mat_fd_unresolved_2026-10-06.rds`

R matrix, 90 sites (rows, SiteID) by 221 taxa (columns, TAXANAME), of mean biomass over the processed replicates (milligrams dry mass per square meter), restricted to the taxa with a trait record and a nonzero total; the matrix passed to FD::dbFD. Copy taken before R/00c first ran.

- Role: input
- Read by: R/00c_resolve_parent_rows.R
- Source: Built upstream from the community table
- One row is: one site; one column per taxon
- R object: table; rows: 90; columns: 221

| Column | Description |
|---|---|
| `Acentrella`, `Acerpenna`, `Acroneuria`, `Aeshnidae`, `Agabus`, `Agnetina` and 215 more (pattern `.*`) | mean biomass over the processed replicates (milligrams dry mass per square meter) of the taxon at the site |

### `den_mat_fd_unresolved_2026-10-06.rds`

R matrix, 90 sites (rows, SiteID) by 221 taxa (columns, TAXANAME), of mean density over the processed replicates (individuals per square meter), restricted to the taxa with a trait record and a nonzero total; the matrix passed to FD::dbFD. Copy taken before R/00c first ran.

- Role: input
- Read by: R/00c_resolve_parent_rows.R
- Source: Built upstream from the community table
- One row is: one site; one column per taxon
- R object: table; rows: 90; columns: 221

| Column | Description |
|---|---|
| `Acentrella`, `Acerpenna`, `Acroneuria`, `Aeshnidae`, `Agabus`, `Agnetina` and 215 more (pattern `.*`) | mean density over the processed replicates (individuals per square meter) of the taxon at the site |

### `hand_drawn_figures/Figure_1.pdf`

Figure 1, the conceptual framework, which has no pipeline source.

- Role: input
- Read by: R/18_manuscript_assets.R (copied into output/manuscript/figures)
- Source: Drawn by the author
- One row is: figure

### `hand_drawn_figures/Figure_1.png`

Figure 1, the conceptual framework, which has no pipeline source.

- Role: input
- Read by: R/18_manuscript_assets.R (copied into output/manuscript/figures)
- Source: Drawn by the author
- One row is: figure

### `hand_drawn_figures/Figure_1.svg`

Figure 1, the conceptual framework, which has no pipeline source.

- Role: input
- Read by: R/18_manuscript_assets.R (copied into output/manuscript/figures)
- Source: Drawn by the author
- One row is: figure

### `site_coordinates.csv`

Location of each stream.

- Role: input
- Read by: R/19_study_area_map.R (Figure 2)
- Source: The three USGS data releases; subset of the project's site file to the 90 Chapter 2 sites
- One row is: one site
- Rows: 90; columns: 5

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Latitude` | Latitude, decimal degrees, WGS 84 |
| `Longitude` | Longitude, decimal degrees, WGS 84 |
| `Region` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `Code` | Two-letter region code: VR, PD or CP |

### `trait_blocks.rds`

R list describing the trait blocks: primary_groups (the 11 CONUS column names), group_labels (the 11 short group codes), trait_cols (the 41 state columns), groups_vec (group index of each state), fuzzy_idx, group_sizes (states per group), gowdis_w (the block-balanced Gower weight of each state, 1 divided by the group's number of states, so every trait group weighs the same in the distance), and state_lookup (state column, group, label).

- Role: input
- Read by: R/00c, R/00a, R/00b, R/00_config.R (IN$trait_blocks: steps 13, 16, 17, 20)
- Source: Step2b
- One row is: not tabular
- R object: list; elements: 8

| Element | Description |
|---|---|
| `primary_groups` | see the paragraph above |
| `group_labels` | see the paragraph above |
| `trait_cols` | see the paragraph above |
| `groups_vec` | see the paragraph above |
| `fuzzy_idx` | see the paragraph above |
| `group_sizes` | see the paragraph above |
| `gowdis_w` | see the paragraph above |
| `state_lookup` | see the paragraph above |

### `traits_binary_proxy_main.rds`

Binary companion of the fuzzy matrix (same taxa, same columns). A trait column is retained in the Gower distance when it is nonzero in either coding, which is why FD::dbFD needs both files to reproduce the published metrics.

- Role: input
- Read by: R/00c, R/00a, R/00b, R/00_config.R (IN$traits_bin: steps 13, 16, 17)
- Source: As traits_fuzzy_proxy_main.rds, with each affinity above zero set to 1
- One row is: one taxon
- R object: table; rows: 576; columns: 43

| Column | Description |
|---|---|
| `Genus` | As in traits_fuzzy_proxy_main.rds |
| `Flying_Strong`, `Flying_Weak` | Adult flying strength: strong or weak; 1 where the fuzzy affinity is above zero, else 0 |
| `Emerge_Fall`, `Emerge_Spring`, `Emerge_Summer`, `Emerge_Winter` | Emergence season: fall, spring, summer or winter; 1 where the fuzzy affinity is above zero, else 0 |
| `Synchro_Poorly`, `Synchro_Well` | Synchronization of emergence: poorly or well synchronized; 1 where the fuzzy affinity is above zero, else 0 |
| `Feed_CF`, `Feed_CG`, `Feed_HB`, `Feed_PA`, `Feed_PR`, `Feed_SH` | Feeding style: CF collector-filterer, CG collector-gatherer, HB herbivore, PA parasite, PR predator, SH shredder; 1 where the fuzzy affinity is above zero, else 0 |
| `Disperse_High`, `Disperse_Low` | Female dispersal: high or low; 1 where the fuzzy affinity is above zero, else 0 |
| `Habit_Burrower`, `Habit_Climber`, `Habit_Clinger`, `Habit_Crawler`, `Habit_Planktonic`, `Habit_Skater` and 2 more (pattern `^Habit_(Burrower|Climber|Clinger|Crawler|Planktonic|Skater|Sprawler|Swimmer)`) | Habit; 1 where the fuzzy affinity is above zero, else 0 |
| `Size_Large`, `Size_Medium`, `Size_Small` | Maximum body size class; 1 where the fuzzy affinity is above zero, else 0 |
| `Respiration_Gills`, `Respiration_Plastron`, `Respiration_Tegument` | Respiration mode; Plastron covers plastron and spiracle (Twardochleb's Plastron_spiracle) and is never merged with gills; 1 where the fuzzy affinity is above zero, else 0 |
| `Rheophily_Both`, `Rheophily_Depositional`, `Rheophily_Erosional` | Rheophily: depositional, erosional or both; 1 where the fuzzy affinity is above zero, else 0 |
| `Thermal_ColdEury`, `Thermal_ColdSteno`, `Thermal_CoolEury`, `Thermal_HotEury`, `Thermal_WarmEury` | Thermal preference: cold stenothermal, cold-cool eurythermal, cool-warm eurythermal, warm eurythermal, hot eurythermal; 1 where the fuzzy affinity is above zero, else 0 |
| `Voltinism_Multivoltine`, `Voltinism_Semivoltine`, `Voltinism_Univoltine` | Voltinism; 1 where the fuzzy affinity is above zero, else 0 |
| `TAXANAME` | Taxon name at the finest rank identified (genus where possible, otherwise family, subfamily, tribe or a higher rank); the key that joins the community table, the matrices and the trait matrix |

### `traits_fuzzy_proxy_main.rds`

Taxon by trait-state affinity matrix, 576 taxa by 41 states in 11 trait groups, built by Data/Data/Step2b_BuildTraitMatrices.R (not shipped; it runs on the full CONUS database) from the CONUS genus trait database of Twardochleb et al. (2021). Affinities within a group sum to 1 where the taxon has a record for that group and are NA where it has none. Taxa identified above genus carry the mean affinity of their closest kin in the database (tribe, then subfamily, then family, then order; Trait_Data/closest_kin_itis.csv), every kin genus with a record counting and a state a genus lacks counting as zero. Seven high-biomass taxa are targeted proxies (Trait_Data/targeted_high_biomass_proxy_*). Codes stand as published; no state is derived from another.

- Role: input
- Read by: R/00c, R/00a, R/00b, R/00_config.R (IN$traits: steps 09, 09b, 11, 15, 15b, 20; IN$traits_fuz: steps 13, 16, 17)
- Source: Twardochleb et al. (2021), CONUS_Genus_Trait_Affinities.csv, through Step2b
- One row is: one taxon
- R object: table; rows: 576; columns: 43

| Column | Description |
|---|---|
| `Genus` | Twardochleb genus, or the community taxon name for a meld or a proxy row |
| `Flying_Strong`, `Flying_Weak` | Adult flying strength: strong or weak; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Emerge_Fall`, `Emerge_Spring`, `Emerge_Summer`, `Emerge_Winter` | Emergence season: fall, spring, summer or winter; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Synchro_Poorly`, `Synchro_Well` | Synchronization of emergence: poorly or well synchronized; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Feed_CF`, `Feed_CG`, `Feed_HB`, `Feed_PA`, `Feed_PR`, `Feed_SH` | Feeding style: CF collector-filterer, CG collector-gatherer, HB herbivore, PA parasite, PR predator, SH shredder; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Disperse_High`, `Disperse_Low` | Female dispersal: high or low; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Habit_Burrower`, `Habit_Climber`, `Habit_Clinger`, `Habit_Crawler`, `Habit_Planktonic`, `Habit_Skater` and 2 more (pattern `^Habit_(Burrower|Climber|Clinger|Crawler|Planktonic|Skater|Sprawler|Swimmer)`) | Habit; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Size_Large`, `Size_Medium`, `Size_Small` | Maximum body size class; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Respiration_Gills`, `Respiration_Plastron`, `Respiration_Tegument` | Respiration mode; Plastron covers plastron and spiracle (Twardochleb's Plastron_spiracle) and is never merged with gills; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Rheophily_Both`, `Rheophily_Depositional`, `Rheophily_Erosional` | Rheophily: depositional, erosional or both; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Thermal_ColdEury`, `Thermal_ColdSteno`, `Thermal_CoolEury`, `Thermal_HotEury`, `Thermal_WarmEury` | Thermal preference: cold stenothermal, cold-cool eurythermal, cool-warm eurythermal, warm eurythermal, hot eurythermal; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `Voltinism_Multivoltine`, `Voltinism_Semivoltine`, `Voltinism_Univoltine` | Voltinism; fuzzy affinity, 0 to 1, summing to 1 within the group |
| `TAXANAME` | Taxon name at the finest rank identified (genus where possible, otherwise family, subfamily, tribe or a higher rank); the key that joins the community table, the matrices and the trait matrix |

## Derived by the pipeline's data steps

### `Regional Analysis/tables/H1/H1_parent_rows_provenance.csv`

The provenance gate (every recomputed column against the unresolved table) and the resolution summary.

- Role: derived; written by R/00c_resolve_parent_rows.R
- Read by: run_all.R (banner), R/11
- Source: R/00c
- One row is: one check
- Rows: 39; columns: 8

| Column | Description |
|---|---|
| `column` | Column checked, or summary quantity |
| `reference` | What it was compared with |
| `max_abs_diff` | Largest absolute difference, or the summary value |
| `n_na` | Missing values |
| `tolerance` | Tolerance applied |
| `reproduced` | TRUE where within tolerance |
| `mode` | gate or resolution summary |
| `stamp` | Date of the unresolved copies |

### `Regional Analysis/tables/H1/H1_parent_rows_resolution.csv`

One row per parent row resolved at a site: what moved to whom.

- Role: derived; written by R/00c_resolve_parent_rows.R
- Read by: R/11, R/18 (output/manuscript/repository_tables/parent_rows_resolved.csv)
- Source: R/00c
- One row is: one parent by site
- Rows: 301; columns: 9

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Typology` | Study region coded as the project's typology number: 1 = Valley and Ridge (Shenandoah), 2 = Coastal Plain (Delmarva), 3 = Piedmont (southeastern Pennsylvania and Maryland) |
| `parent` | The parent row's taxon name |
| `rank` | Its rank: family, subfamily, tribe, superfamily or order_or_higher |
| `n_children` | Children present at the site |
| `children` | Their names, separated by semicolons |
| `abundance_moved` | Corrected abundance distributed |
| `density_moved` | Density distributed |
| `biomass_moved` | Biomass distributed |

### `Regional Analysis/tables/H1/H1_redundancy_geometry_provenance.csv`

Checks that the binary geometry reproduces the upstream redundancy columns before the fuzzy geometry replaces them.

- Role: derived; written by R/00b_redundancy_geometry.R
- Read by: no step
- Source: R/00b
- One row is: one check
- Rows: 9; columns: 2

| Column | Description |
|---|---|
| `check` | What was checked |
| `value` | Result |

### `Regional Analysis/tables/H1/H1_site_summary_macro_metrics.csv`

The site table after parent-row resolution: R/00c's recomputed metrics with the unresolved table's other columns carried.

- Role: derived; written by R/00c_resolve_parent_rows.R
- Read by: R/00a_trait_columns.R
- Source: R/00c from the resolved community table
- One row is: one site
- Rows: 90; columns: 168

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Typology` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `FRic_den_binPool`, `FEve_den_binPool`, `FDiv_den_binPool`, `FDis_den_binPool`, `RaoQ_den_binPool`, `FRic_bio_binPool` and 33 more (pattern `^(FRic|FEve|FDiv|FDis|RaoQ)_(den|bio)_(binPool|fuzPool|binReg|fuzReg)$`) | Legacy variant of the metric: binary (bin) or fuzzy (fuz) coding, in the pooled (Pool) or regional (Reg) trait space; dropped by R/00a except FRic_den_fuzPool |
| `FRic_den_fuzPool` | Functional richness computed in the trait space of the whole taxon pool (all three regions) under the fuzzy coding; read by R/00b for space-filling redundancy |
| `FRic_den` | Functional richness: convex hull volume of the site's trait set in the ten-axis PCoA of the fuzzy group-weighted Gower distance (FD::dbFD, Cailliez correction); presence-based, one value per site |
| `FEve_den`, `FDiv_den`, `FDis_den`, `RaoQ_den`, `FEve_bio`, `FDiv_bio` and 2 more (pattern `^(FEve|FDiv|FDis|RaoQ)_(den|bio)$`) | Functional evenness, divergence, dispersion or Rao's quadratic entropy (FD::dbFD on the fuzzy proxy trait matrix with block-balanced Gower weights, Cailliez correction, ten axes), weighted by mean density (den) or mean biomass (bio) over the processed replicates; Villeger et al. 2008, Laliberte and Legendre 2010, Botta-Dukat 2005 |
| `Simpson_den` | Simpson diversity 1 - sum(p_i^2) on mean density, over the taxa with a trait record |
| `Simpson_bio` | Simpson diversity 1 - sum(p_i^2) on mean biomass, over the taxa with a trait record |
| `Richness_FD` | Number of taxa with a trait record and a nonzero total at the site (the trait set) |
| `RaoQ_gower_den` | Rao's quadratic entropy on the binary unweighted Gower distance, density weights; the geometry the upstream table used for redundancy, replaced by R/00b's fuzzy geometry |
| `RaoQ_gower_bio` | As RaoQ_gower_den, biomass weights |
| `FR_RD_den` | Overlap redundancy weighted by density, 1 - Q/D, with Q Rao's quadratic entropy on the fuzzy group-weighted Gower distance and D Simpson diversity (de Bello et al. 2007; Ricotta et al. 2016) |
| `FR_RD_bio` | Overlap redundancy weighted by biomass, 1 - Q/D |
| `FR_SF_den` | Space-filling redundancy: Richness_FD divided by FRic_den_fuzPool (Mouillot et al. 2014) |
| `FR_SF_den_log10` | log10 of FR_SF_den; the screened form |
| `FR_DB_RD` | Redundancy divergence: FR_RD_bio minus FR_RD_den |
| `Richness_alpha` | Taxonomic richness: number of taxa present at the site, at the finest rank identified |
| `Shannon` | Shannon diversity H' on mean density over the processed replicates, natural logarithm, over every taxon |
| `Total_Density` | Total density, individuals per square meter: the median over the site's processed replicate samples of each replicate's total (R/00d, from Step2_SiteReplicateMedians.csv; the lab method M3). On the unresolved tables, before R/00d, the mean over the processed replicates |
| `Total_Density_mean` | Total density as the mean over the processed replicates (the sum over taxa of the mean taxon density); shipped beside the median so both site-total bases are reported (Table S1). Written by R/00d |
| `Total_Biomass` | Total biomass, milligrams dry mass per square meter: the median over the processed replicates of each replicate's total (R/00d; M3). On the unresolved tables, the mean over the processed replicates |
| `Total_Biomass_mean` | Total biomass as the mean over the processed replicates; written by R/00d |
| `BIBI_score` | Chessie B-IBI score of the site as the upstream site table carried it; the pipeline reads the index from H1_BIBI_site_with_classes_MEDIAN.csv instead, so this column is carried, never used |
| `BIBI_Class` | Narrative rating carried from the upstream site table; not used by the pipeline (see BIBI_score) |
| `Disperse_High.dispersal_CWM_den`, `Disperse_Low.dispersal_CWM_den`, `Emerge_Spring_CWM_den`, `Emerge_Summer_CWM_den`, `Emerge_Winter_CWM_den`, `Feed_CF_Collector-filterer_CWM_den` and 70 more (pattern `^(Flying|Emerge|Synchro|Feed|Disperse|Habit|Size|Respiration|Rheophily|Thermal|Voltinism)_.*_CWM_(den|bio)$`) | Legacy community-weighted mean from the binary trait file (a taxon with several states counted in full in each); dropped by R/00a |
| `CWM_Response_Density`, `CWM_Response_Biomass`, `CWM_Effect_Density`, `CWM_Effect_Biomass`, `CWM_Linking_Density`, `CWM_Linking_Biomass` and 6 more (pattern `^(CWM|REL)_(Effect|Linking|Response)_(Density|Biomass|den|bio)$`) | Legacy aggregate of the binary trait file by trait role (effect, linking, response); dropped by R/00a |
| `HBI_taxa_with_TV` | Number of taxa at the site with a tolerance value |
| `HBI_taxa_total` | Number of taxa at the site |
| `HBI_rel_abund` | Hilsenhoff Biotic Index weighted by relative abundance (the conventional form; family tolerance values from the Chessie B-IBI, Buchanan et al. 2011); the median over the processed replicates of each replicate's index after R/00d, the index on the mean taxon values before it |
| `HBI_rel_abund_mean` | HBI_rel_abund computed on the mean taxon values (the version before R/00d); written by R/00d |
| `HBI_den` | Hilsenhoff Biotic Index weighted by density; identical to HBI_rel_abund within a site and not screened |
| `HBI_den_mean` | HBI_den computed on the mean taxon values; written by R/00d |
| `HBI_bio` | Hilsenhoff Biotic Index weighted by biomass (median over the processed replicates after R/00d) |
| `HBI_bio_mean` | HBI_bio computed on the mean taxon values; written by R/00d |
| `HBI_coverage_pct` | Percent of abundance in taxa with a tolerance value |
| `HBI_density_coverage_pct` | Percent of density in taxa with a tolerance value |
| `HBI_biomass_coverage_pct` | Percent of biomass in taxa with a tolerance value |

### `Regional Analysis/tables/H1/H1_site_summary_macro_metrics_fuzzyRD.csv`

The site table run_all.R reads: the redundancy family (RaoQ_gower, FR_RD, FR_DB_RD, FRic_den, FR_SF) recomputed in the fuzzy group-weighted Gower geometry, every other column as R/00a left it.

- Role: derived; written by R/00b_redundancy_geometry.R
- Read by: R/00_config.R (IN$site_metrics under CH2_RD_GEOMETRY=fuzzy): every step
- Source: R/00b on R/00a's output
- One row is: one site
- Rows: 90; columns: 123

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Typology` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `FRic_den_fuzPool` | Functional richness computed in the trait space of the whole taxon pool (all three regions) under the fuzzy coding; read by R/00b for space-filling redundancy |
| `FRic_den` | Functional richness: convex hull volume of the site's trait set in the ten-axis PCoA of the fuzzy group-weighted Gower distance (FD::dbFD, Cailliez correction); presence-based, one value per site |
| `FEve_den`, `FDiv_den`, `FDis_den`, `RaoQ_den`, `FEve_bio`, `FDiv_bio` and 2 more (pattern `^(FEve|FDiv|FDis|RaoQ)_(den|bio)$`) | Functional evenness, divergence, dispersion or Rao's quadratic entropy (FD::dbFD on the fuzzy proxy trait matrix with block-balanced Gower weights, Cailliez correction, ten axes), weighted by mean density (den) or mean biomass (bio) over the processed replicates; Villeger et al. 2008, Laliberte and Legendre 2010, Botta-Dukat 2005 |
| `Simpson_den` | Simpson diversity 1 - sum(p_i^2) on mean density, over the taxa with a trait record |
| `Simpson_bio` | Simpson diversity 1 - sum(p_i^2) on mean biomass, over the taxa with a trait record |
| `Richness_FD` | Number of taxa with a trait record and a nonzero total at the site (the trait set) |
| `RaoQ_gower_den` | Rao's quadratic entropy on the binary unweighted Gower distance, density weights; the geometry the upstream table used for redundancy, replaced by R/00b's fuzzy geometry |
| `RaoQ_gower_bio` | As RaoQ_gower_den, biomass weights |
| `FR_RD_den` | Overlap redundancy weighted by density, 1 - Q/D, with Q Rao's quadratic entropy on the fuzzy group-weighted Gower distance and D Simpson diversity (de Bello et al. 2007; Ricotta et al. 2016) |
| `FR_RD_bio` | Overlap redundancy weighted by biomass, 1 - Q/D |
| `FR_SF_den` | Space-filling redundancy: Richness_FD divided by FRic_den_fuzPool (Mouillot et al. 2014) |
| `FR_SF_den_log10` | log10 of FR_SF_den; the screened form |
| `FR_DB_RD` | Redundancy divergence: FR_RD_bio minus FR_RD_den |
| `Richness_alpha` | Taxonomic richness: number of taxa present at the site, at the finest rank identified |
| `Shannon` | Shannon diversity H' on mean density over the processed replicates, natural logarithm, over every taxon |
| `Total_Density` | Total density, individuals per square meter: the median over the site's processed replicate samples of each replicate's total (R/00d, from Step2_SiteReplicateMedians.csv; the lab method M3). On the unresolved tables, before R/00d, the mean over the processed replicates |
| `Total_Density_mean` | Total density as the mean over the processed replicates (the sum over taxa of the mean taxon density); shipped beside the median so both site-total bases are reported (Table S1). Written by R/00d |
| `Total_Biomass` | Total biomass, milligrams dry mass per square meter: the median over the processed replicates of each replicate's total (R/00d; M3). On the unresolved tables, the mean over the processed replicates |
| `Total_Biomass_mean` | Total biomass as the mean over the processed replicates; written by R/00d |
| `BIBI_score` | Chessie B-IBI score of the site as the upstream site table carried it; the pipeline reads the index from H1_BIBI_site_with_classes_MEDIAN.csv instead, so this column is carried, never used |
| `BIBI_Class` | Narrative rating carried from the upstream site table; not used by the pipeline (see BIBI_score) |
| `HBI_taxa_with_TV` | Number of taxa at the site with a tolerance value |
| `HBI_taxa_total` | Number of taxa at the site |
| `HBI_rel_abund` | Hilsenhoff Biotic Index weighted by relative abundance (the conventional form; family tolerance values from the Chessie B-IBI, Buchanan et al. 2011); the median over the processed replicates of each replicate's index after R/00d, the index on the mean taxon values before it |
| `HBI_rel_abund_mean` | HBI_rel_abund computed on the mean taxon values (the version before R/00d); written by R/00d |
| `HBI_den` | Hilsenhoff Biotic Index weighted by density; identical to HBI_rel_abund within a site and not screened |
| `HBI_den_mean` | HBI_den computed on the mean taxon values; written by R/00d |
| `HBI_bio` | Hilsenhoff Biotic Index weighted by biomass (median over the processed replicates after R/00d) |
| `HBI_bio_mean` | HBI_bio computed on the mean taxon values; written by R/00d |
| `HBI_coverage_pct` | Percent of abundance in taxa with a tolerance value |
| `HBI_density_coverage_pct` | Percent of density in taxa with a tolerance value |
| `HBI_biomass_coverage_pct` | Percent of biomass in taxa with a tolerance value |
| `Flying_Strong_CWM_den`, `Flying_Weak_CWM_den`, `Flying_Strong_CWM_bio`, `Flying_Weak_CWM_bio` | Adult flying strength: strong or weak. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Emerge_Fall_CWM_den`, `Emerge_Spring_CWM_den`, `Emerge_Summer_CWM_den`, `Emerge_Winter_CWM_den`, `Emerge_Fall_CWM_bio`, `Emerge_Spring_CWM_bio` and 2 more (pattern `^Emerge_(Fall|Spring|Summer|Winter)_CWM_(den|bio)$`) | Emergence season: fall, spring, summer or winter. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Synchro_Poorly_CWM_den`, `Synchro_Well_CWM_den`, `Synchro_Poorly_CWM_bio`, `Synchro_Well_CWM_bio` | Synchronization of emergence: poorly or well synchronized. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Feed_CF_CWM_den`, `Feed_CG_CWM_den`, `Feed_HB_CWM_den`, `Feed_PA_CWM_den`, `Feed_PR_CWM_den`, `Feed_SH_CWM_den` and 6 more (pattern `^Feed_(CF|CG|HB|PA|PR|SH)_CWM_(den|bio)$`) | Feeding style: CF collector-filterer, CG collector-gatherer, HB herbivore, PA parasite, PR predator, SH shredder. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Disperse_High_CWM_den`, `Disperse_Low_CWM_den`, `Disperse_High_CWM_bio`, `Disperse_Low_CWM_bio` | Female dispersal: high or low. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Habit_Burrower_CWM_den`, `Habit_Climber_CWM_den`, `Habit_Clinger_CWM_den`, `Habit_Crawler_CWM_den`, `Habit_Planktonic_CWM_den`, `Habit_Skater_CWM_den` and 10 more (pattern `^Habit_(Burrower|Climber|Clinger|Crawler|Planktonic|Skater|Sprawler|Swimmer)_CWM_(den|bio)$`) | Habit. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Size_Large_CWM_den`, `Size_Medium_CWM_den`, `Size_Small_CWM_den`, `Size_Large_CWM_bio`, `Size_Medium_CWM_bio`, `Size_Small_CWM_bio` | Maximum body size class. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Respiration_Gills_CWM_den`, `Respiration_Plastron_CWM_den`, `Respiration_Tegument_CWM_den`, `Respiration_Gills_CWM_bio`, `Respiration_Plastron_CWM_bio`, `Respiration_Tegument_CWM_bio` | Respiration mode; Plastron covers plastron and spiracle (Twardochleb's Plastron_spiracle) and is never merged with gills. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Rheophily_Both_CWM_den`, `Rheophily_Depositional_CWM_den`, `Rheophily_Erosional_CWM_den`, `Rheophily_Both_CWM_bio`, `Rheophily_Depositional_CWM_bio`, `Rheophily_Erosional_CWM_bio` | Rheophily: depositional, erosional or both. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Thermal_ColdEury_CWM_den`, `Thermal_ColdSteno_CWM_den`, `Thermal_CoolEury_CWM_den`, `Thermal_HotEury_CWM_den`, `Thermal_WarmEury_CWM_den`, `Thermal_ColdEury_CWM_bio` and 4 more (pattern `^Thermal_(ColdEury|ColdSteno|CoolEury|HotEury|WarmEury)_CWM_(den|bio)$`) | Thermal preference: cold stenothermal, cold-cool eurythermal, cool-warm eurythermal, warm eurythermal, hot eurythermal. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Voltinism_Multivoltine_CWM_den`, `Voltinism_Semivoltine_CWM_den`, `Voltinism_Univoltine_CWM_den`, `Voltinism_Multivoltine_CWM_bio`, `Voltinism_Semivoltine_CWM_bio`, `Voltinism_Univoltine_CWM_bio` | Voltinism. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |

### `Regional Analysis/tables/H1/H1_site_summary_macro_metrics_twardochleb.csv`

The resolved site table with the functional diversity and community-weighted mean columns rebuilt from the Twardochleb fuzzy matrix and the legacy columns dropped.

- Role: derived; written by R/00a_trait_columns.R
- Read by: R/00b_redundancy_geometry.R
- Source: R/00a on the resolved data
- One row is: one site
- Rows: 90; columns: 118

| Column | Description |
|---|---|
| `SiteID` | Stream code, four letters, as in the U.S. Geological Survey data releases |
| `Typology` | Study region by name: Valley & Ridge, Piedmont or Coastal Plain |
| `FRic_den_fuzPool` | Functional richness computed in the trait space of the whole taxon pool (all three regions) under the fuzzy coding; read by R/00b for space-filling redundancy |
| `FRic_den` | Functional richness: convex hull volume of the site's trait set in the ten-axis PCoA of the fuzzy group-weighted Gower distance (FD::dbFD, Cailliez correction); presence-based, one value per site |
| `FEve_den`, `FDiv_den`, `FDis_den`, `RaoQ_den`, `FEve_bio`, `FDiv_bio` and 2 more (pattern `^(FEve|FDiv|FDis|RaoQ)_(den|bio)$`) | Functional evenness, divergence, dispersion or Rao's quadratic entropy (FD::dbFD on the fuzzy proxy trait matrix with block-balanced Gower weights, Cailliez correction, ten axes), weighted by mean density (den) or mean biomass (bio) over the processed replicates; Villeger et al. 2008, Laliberte and Legendre 2010, Botta-Dukat 2005 |
| `Simpson_den` | Simpson diversity 1 - sum(p_i^2) on mean density, over the taxa with a trait record |
| `Simpson_bio` | Simpson diversity 1 - sum(p_i^2) on mean biomass, over the taxa with a trait record |
| `Richness_FD` | Number of taxa with a trait record and a nonzero total at the site (the trait set) |
| `RaoQ_gower_den` | Rao's quadratic entropy on the binary unweighted Gower distance, density weights; the geometry the upstream table used for redundancy, replaced by R/00b's fuzzy geometry |
| `RaoQ_gower_bio` | As RaoQ_gower_den, biomass weights |
| `FR_RD_den` | Overlap redundancy weighted by density, 1 - Q/D, with Q Rao's quadratic entropy on the fuzzy group-weighted Gower distance and D Simpson diversity (de Bello et al. 2007; Ricotta et al. 2016) |
| `FR_RD_bio` | Overlap redundancy weighted by biomass, 1 - Q/D |
| `FR_SF_den` | Space-filling redundancy: Richness_FD divided by FRic_den_fuzPool (Mouillot et al. 2014) |
| `FR_SF_den_log10` | log10 of FR_SF_den; the screened form |
| `FR_DB_RD` | Redundancy divergence: FR_RD_bio minus FR_RD_den |
| `Richness_alpha` | Taxonomic richness: number of taxa present at the site, at the finest rank identified |
| `Shannon` | Shannon diversity H' on mean density over the processed replicates, natural logarithm, over every taxon |
| `Total_Density` | Total density, individuals per square meter: the median over the site's processed replicate samples of each replicate's total (R/00d, from Step2_SiteReplicateMedians.csv; the lab method M3). On the unresolved tables, before R/00d, the mean over the processed replicates |
| `Total_Biomass` | Total biomass, milligrams dry mass per square meter: the median over the processed replicates of each replicate's total (R/00d; M3). On the unresolved tables, the mean over the processed replicates |
| `BIBI_score` | Chessie B-IBI score of the site as the upstream site table carried it; the pipeline reads the index from H1_BIBI_site_with_classes_MEDIAN.csv instead, so this column is carried, never used |
| `BIBI_Class` | Narrative rating carried from the upstream site table; not used by the pipeline (see BIBI_score) |
| `HBI_taxa_with_TV` | Number of taxa at the site with a tolerance value |
| `HBI_taxa_total` | Number of taxa at the site |
| `HBI_rel_abund` | Hilsenhoff Biotic Index weighted by relative abundance (the conventional form; family tolerance values from the Chessie B-IBI, Buchanan et al. 2011); the median over the processed replicates of each replicate's index after R/00d, the index on the mean taxon values before it |
| `HBI_den` | Hilsenhoff Biotic Index weighted by density; identical to HBI_rel_abund within a site and not screened |
| `HBI_bio` | Hilsenhoff Biotic Index weighted by biomass (median over the processed replicates after R/00d) |
| `HBI_coverage_pct` | Percent of abundance in taxa with a tolerance value |
| `HBI_density_coverage_pct` | Percent of density in taxa with a tolerance value |
| `HBI_biomass_coverage_pct` | Percent of biomass in taxa with a tolerance value |
| `Flying_Strong_CWM_den`, `Flying_Weak_CWM_den`, `Flying_Strong_CWM_bio`, `Flying_Weak_CWM_bio` | Adult flying strength: strong or weak. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Emerge_Fall_CWM_den`, `Emerge_Spring_CWM_den`, `Emerge_Summer_CWM_den`, `Emerge_Winter_CWM_den`, `Emerge_Fall_CWM_bio`, `Emerge_Spring_CWM_bio` and 2 more (pattern `^Emerge_(Fall|Spring|Summer|Winter)_CWM_(den|bio)$`) | Emergence season: fall, spring, summer or winter. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Synchro_Poorly_CWM_den`, `Synchro_Well_CWM_den`, `Synchro_Poorly_CWM_bio`, `Synchro_Well_CWM_bio` | Synchronization of emergence: poorly or well synchronized. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Feed_CF_CWM_den`, `Feed_CG_CWM_den`, `Feed_HB_CWM_den`, `Feed_PA_CWM_den`, `Feed_PR_CWM_den`, `Feed_SH_CWM_den` and 6 more (pattern `^Feed_(CF|CG|HB|PA|PR|SH)_CWM_(den|bio)$`) | Feeding style: CF collector-filterer, CG collector-gatherer, HB herbivore, PA parasite, PR predator, SH shredder. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Disperse_High_CWM_den`, `Disperse_Low_CWM_den`, `Disperse_High_CWM_bio`, `Disperse_Low_CWM_bio` | Female dispersal: high or low. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Habit_Burrower_CWM_den`, `Habit_Climber_CWM_den`, `Habit_Clinger_CWM_den`, `Habit_Crawler_CWM_den`, `Habit_Planktonic_CWM_den`, `Habit_Skater_CWM_den` and 10 more (pattern `^Habit_(Burrower|Climber|Clinger|Crawler|Planktonic|Skater|Sprawler|Swimmer)_CWM_(den|bio)$`) | Habit. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Size_Large_CWM_den`, `Size_Medium_CWM_den`, `Size_Small_CWM_den`, `Size_Large_CWM_bio`, `Size_Medium_CWM_bio`, `Size_Small_CWM_bio` | Maximum body size class. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Respiration_Gills_CWM_den`, `Respiration_Plastron_CWM_den`, `Respiration_Tegument_CWM_den`, `Respiration_Gills_CWM_bio`, `Respiration_Plastron_CWM_bio`, `Respiration_Tegument_CWM_bio` | Respiration mode; Plastron covers plastron and spiracle (Twardochleb's Plastron_spiracle) and is never merged with gills. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Rheophily_Both_CWM_den`, `Rheophily_Depositional_CWM_den`, `Rheophily_Erosional_CWM_den`, `Rheophily_Both_CWM_bio`, `Rheophily_Depositional_CWM_bio`, `Rheophily_Erosional_CWM_bio` | Rheophily: depositional, erosional or both. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Thermal_ColdEury_CWM_den`, `Thermal_ColdSteno_CWM_den`, `Thermal_CoolEury_CWM_den`, `Thermal_HotEury_CWM_den`, `Thermal_WarmEury_CWM_den`, `Thermal_ColdEury_CWM_bio` and 4 more (pattern `^Thermal_(ColdEury|ColdSteno|CoolEury|HotEury|WarmEury)_CWM_(den|bio)$`) | Thermal preference: cold stenothermal, cold-cool eurythermal, cool-warm eurythermal, warm eurythermal, hot eurythermal. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |
| `Voltinism_Multivoltine_CWM_den`, `Voltinism_Semivoltine_CWM_den`, `Voltinism_Univoltine_CWM_den`, `Voltinism_Multivoltine_CWM_bio`, `Voltinism_Semivoltine_CWM_bio`, `Voltinism_Univoltine_CWM_bio` | Voltinism. Community-weighted mean affinity per site, weighted by density (den) or biomass (bio), over the taxa with a record for the group |

### `Regional Analysis/tables/H1/H1_trait_columns_provenance.csv`

For each rebuilt functional diversity column, how closely R/00a's recomputation matched the upstream value.

- Role: derived; written by R/00a_trait_columns.R
- Read by: no step
- Source: R/00a
- One row is: one column
- Rows: 15; columns: 6

| Column | Description |
|---|---|
| `column` | Column |
| `max_abs_diff_from_upstream` | Largest absolute difference |
| `spearman_with_upstream` | Spearman correlation with the upstream column |
| `n_nonfinite` | Non-finite values |
| `mode` | rebuild or verify |
| `note` | Note |

### `SiteTaxa_Density_Biomass_Transformed.csv`

The community table after parent rows are resolved by DPAC (distribute parent among children; Cuffney et al. 2007): at each site, a family, subfamily, tribe or higher-rank row whose children occur at that site gives its abundance, density and biomass to those children in proportion to their density and leaves the site. Same columns as the unresolved copy.

- Role: derived; written by R/00c_resolve_parent_rows.R
- Read by: R/10, R/11, R/12 (through the site table and matrices, every step)
- Source: Written by R/00c from the unresolved copy
- One row is: one site by taxon
- Rows: 3740; columns: 26

| Column | Description |
|---|---|
| `Typology`, `SiteID`, `TAXANAME`, `Family`, `Mean_Corrected_Abundance`, `Median_Corrected_Abundance` and 20 more (pattern `.*`) | As in the unresolved copy |

### `bio_mat_fd.rds`

R matrix, 90 sites by 215 taxa, of mean biomass over the processed replicates (milligrams dry mass per square meter) after parent rows are resolved.

- Role: derived; written by R/00c_resolve_parent_rows.R
- Read by: R/00a, R/00b, R/00_config.R (every step)
- Source: Written by R/00c from the resolved community table
- One row is: one site; one column per taxon
- R object: table; rows: 90; columns: 215

| Column | Description |
|---|---|
| `Acentrella`, `Acerpenna`, `Acroneuria`, `Aeshnidae`, `Agabus`, `Agnetina` and 209 more (pattern `.*`) | mean biomass over the processed replicates (milligrams dry mass per square meter) of the taxon at the site |

### `den_mat_fd.rds`

R matrix, 90 sites by 215 taxa, of mean density over the processed replicates (individuals per square meter) after parent rows are resolved.

- Role: derived; written by R/00c_resolve_parent_rows.R
- Read by: R/00a, R/00b, R/00_config.R (every step)
- Source: Written by R/00c from the resolved community table
- One row is: one site; one column per taxon
- R object: table; rows: 90; columns: 215

| Column | Description |
|---|---|
| `Acentrella`, `Acerpenna`, `Acroneuria`, `Aeshnidae`, `Agabus`, `Agnetina` and 209 more (pattern `.*`) | mean density over the processed replicates (individuals per square meter) of the taxon at the site |

## Documentation

### `Trait_Data/family_gap_fill_log.csv`

Every trait group filled from closest kin: which taxon, which group, which kin rank, how many kin genera had a record.

- Role: documentation
- Read by: no pipeline step
- Source: Step2b
- One row is: one taxon by trait group
- Rows: 869; columns: 8

| Column | Description |
|---|---|
| `matrix` | main or sensitivity |
| `TAXANAME` | Taxon name at the finest rank identified (genus where possible, otherwise family, subfamily, tribe or a higher rank); the key that joins the community table, the matrices and the trait matrix |
| `row_kind` | Twardochleb genus, meld or proxy |
| `family` | Family of the taxon |
| `kin_rank_used` | Rank at which kin were pooled |
| `trait_group` | Trait group filled |
| `action` | filled from closest kin, or left NA |
| `n_kin_genera_with_record` | Kin genera with a record for the group |

### `Trait_Data/genus_family_resolved.csv`

The one family assigned to each genus for the kin rule, with the families the database listed for it.

- Role: documentation
- Read by: Step2b (upstream); no pipeline step
- Source: ITIS, accessed 22 September 2026
- One row is: one genus
- Rows: 18; columns: 4

| Column | Description |
|---|---|
| `Genus` | Genus |
| `Family` | Family used |
| `families_listed` | Families the database listed, with counts |
| `source` | ITIS serial number and access date |

### `Trait_Data/targeted_high_biomass_proxy_manifest.csv`

The seven high-biomass taxa given a targeted proxy, with the genera that filled each missing trait group.

- Role: documentation
- Read by: no pipeline step
- Source: Step2b and targeted_high_biomass_proxy_rationale.md
- One row is: one taxon
- Rows: 7; columns: 13

| Column | Description |
|---|---|
| `TAXANAME` | Taxon name at the finest rank identified (genus where possible, otherwise family, subfamily, tribe or a higher rank); the key that joins the community table, the matrices and the trait matrix |
| `ProxyTier` | targeted_high_biomass |
| `ProxyMethod` | exact_plus_targeted_fill or targeted_fill |
| `IncludeMain` | TRUE |
| `IncludeSensitivity` | TRUE |
| `Confidence` | Author's confidence in the assignment |
| `ExactSource` | Twardochleb genus matched exactly, where one exists |
| `FillSourceGenera` | Genera whose states filled missing groups |
| `FillTaxaFound` | Of those, the genera found in the database |
| `n_trait_groups_filled` | Trait groups recovered |
| `ProxyStatus` | included_targeted_proxy |
| `Rationale` | Why these genera |
| `GroupSourceAudit` | Per group, exact or fill |

### `Trait_Data/targeted_high_biomass_proxy_rationale.md`

Prose rationale for the targeted proxies.

- Role: documentation
- Read by: no pipeline step
- Source: Author
- One row is: prose

### `Trait_Data/trait_proxy_manifest.csv`

One row per community taxon whose trait row was built by proxy: the method, the source taxon or genera, the number of source genera, the tier, and whether it enters the main and the sensitivity matrix.

- Role: documentation
- Read by: no pipeline step
- Source: Step2b
- One row is: one proxy taxon
- Rows: 100; columns: 20

| Column | Description |
|---|---|
| `TAXANAME` | Taxon name at the finest rank identified (genus where possible, otherwise family, subfamily, tribe or a higher rank); the key that joins the community table, the matrices and the trait matrix |
| `method` | How the row was built |
| `proxy_source` | Taxon whose genera were pooled |
| `n_source_genera` | Number of source genera |
| `tier` | Proxy tier (1 exact, 2 family or subfamily aggregate, 3 order) |
| `included_main` | TRUE where the row enters the main matrix |
| `included_sensitivity` | TRUE where it enters the sensitivity matrix |
| `notes` | Free text |
| `ProxyTier`, `ProxyMethod`, `IncludeMain`, `IncludeSensitivity`, `Confidence`, `ExactSource` and 6 more (pattern `^(ProxyTier|ProxyMethod|IncludeMain|IncludeSensitivity|Confidence|ExactSource|FillSourceGenera|FillTaxaFound|n_trait_groups_filled|ProxyStatus|Rationale|GroupSourceAudit)$`) | Filled only for the targeted high-biomass proxies; see targeted_high_biomass_proxy_manifest.csv |

